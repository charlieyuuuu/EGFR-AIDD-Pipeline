"""
03_qsar_model.py
Description: Build Random Forest QSAR model using Morgan Fingerprints.
"""
nano run_qsar.py

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')

# ==================== 第一步：准备数据 ====================
print("1. 加载数据...")
df = pd.read_csv('CHEMBL203_actives.csv')
print(f"   原始数据: {len(df)} 条")

# 清理数据：去掉没有 SMILES 或 pChEMBL 值的行
df = df.dropna(subset=['smiles', 'pChEMBL'])
df['pChEMBL'] = pd.to_numeric(df['pChEMBL'], errors='coerce')
df = df.dropna(subset=['pChEMBL'])
print(f"   清理后: {len(df)} 条")

# ==================== 第二步：生成分子指纹 ====================
print("2. 生成 Morgan 指纹 (ECFP4)...")
def smiles_to_fingerprint(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    # 生成 2048 位的 Morgan 指纹，半径 2 对应 ECFP4
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return np.array(fp)

# 对每个分子生成指纹，记录哪些失败了
fingerprints = []
labels = []
for idx, row in df.iterrows():
    fp = smiles_to_fingerprint(row['smiles'])
    if fp is not None:
        fingerprints.append(fp)
        labels.append(row['pChEMBL'])

X = np.array(fingerprints)
y = np.array(labels)
print(f"   成功生成指纹: {X.shape[0]} 个分子, 每个 {X.shape[1]} 维特征")

# ==================== 第三步：划分数据并训练模型 ====================
print("3. 划分训练集/测试集并训练随机森林...")
# 80% 训练，20% 测试
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)
print(f"   训练集: {len(X_train)} 条, 测试集: {len(X_test)} 条")

# 随机森林回归模型：用 100 棵树，用所有 CPU 核心
model = RandomForestRegressor(
    n_estimators=100,
    max_depth=20,
    min_samples_leaf=3,
    random_state=42,
    n_jobs=-1
)
model.fit(X_train, y_train)

# ==================== 第四步：评估模型 ====================
print("4. 评估模型性能...")
y_pred_train = model.predict(X_train)
y_pred_test = model.predict(X_test)

r2_train = r2_score(y_train, y_pred_train)
r2_test = r2_score(y_test, y_pred_test)
rmse_test = np.sqrt(mean_squared_error(y_test, y_pred_test))

print(f"   训练集 R²: {r2_train:.3f}")
print(f"   测试集 R²: {r2_test:.3f}")
print(f"   测试集 RMSE: {rmse_test:.3f} (pIC50 单位)")
print("   (R² 越接近 1 越好，0.5-0.7 在 QSAR 里已经算可用)")

# ==================== 第五步：预测你的 top50 ====================
print("5. 预测 top50 对接结果的活性...")
# 假设你已经生成了 top50_docking.csv，包含 chembl_id 和 affinity
# 如果没有，可以从 docking_results.csv 取前 50
try:
    top50 = pd.read_csv('top50_docking.csv')
except FileNotFoundError:
    top50 = pd.read_csv('docking_results.csv').sort_values('affinity').head(50)
    top50.to_csv('top50_docking.csv', index=False)

# 从原始数据里找到这些化合物的 SMILES
# top50 里的 chembl_id 是 CHEMBL 编号，从 df 里匹配
df_lookup = df.set_index('chembl_id')['smiles'].to_dict()

top_fps = []
top_valid = []
for _, row in top50.iterrows():
    cid = row['chembl_id']
    smi = df_lookup.get(cid)
    if smi is None:
        continue
    fp = smiles_to_fingerprint(smi)
    if fp is not None:
        top_fps.append(fp)
        top_valid.append({
            'chembl_id': cid,
            'affinity': row['affinity'],
            'smiles': smi
        })

if top_fps:
    top_X = np.array(top_fps)
    top_pred = model.predict(top_X)
    for i, item in enumerate(top_valid):
        item['predicted_pIC50'] = round(top_pred[i], 2)

    result_df = pd.DataFrame(top_valid)
    result_df = result_df.sort_values('predicted_pIC50', ascending=False)
    result_df.to_csv('top50_qsar_prediction.csv', index=False)
    print(f"   已预测 {len(result_df)} 个化合物，结果保存到 top50_qsar_prediction.csv")
    print(result_df.head(10).to_string(index=False))
else:
    print("   top50 里没有找到可用的 SMILES，请检查 chembl_id 是否匹配")

nohup python run_qsar.py > run_qsar.log 2>&1 &

nano correlation_analysis.py

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 读取 top50 的 QSAR 预测结果
df = pd.read_csv('top50_qsar_prediction.csv')
print(f"化合物数: {len(df)}")

# 对接打分越负越好，pIC50 越高越好
# 为了直观，把 affinity 取反，变成"对接得分"（越高越好）
df['docking_score'] = -df['affinity']

# 计算相关系数
pearson_r, pearson_p = pearsonr(df['docking_score'], df['predicted_pIC50'])
spearman_r, spearman_p = spearmanr(df['docking_score'], df['predicted_pIC50'])

print(f"\nPearson 相关系数: {pearson_r:.3f} (p={pearson_p:.4f})")
print(f"Spearman 相关系数: {spearman_r:.3f} (p={spearman_p:.4f})")
print("\n解释: 相关系数 > 0 表示对接打分越好，QSAR 预测活性也越高")

# 找矛盾分子：对接打分高但 QSAR 预测低，或反过来
df['rank_docking'] = df['docking_score'].rank(ascending=False)
df['rank_qsar'] = df['predicted_pIC50'].rank(ascending=False)
df['rank_diff'] = abs(df['rank_docking'] - df['rank_qsar'])

print("\n排名差异最大的 10 个（对接和 QSAR 分歧最大）:")
print(df.nlargest(10, 'rank_diff')[
    ['chembl_id', 'affinity', 'predicted_pIC50', 'rank_docking', 'rank_qsar']
].to_string(index=False))

# 画散点图
plt.figure(figsize=(8, 6))
plt.scatter(df['docking_score'], df['predicted_pIC50'], alpha=0.7, edgecolors='k')
plt.xlabel('Docking Score (-affinity, kcal/mol)')
plt.ylabel('Predicted pIC50')
plt.title(f'Docking Score vs QSAR Predicted Activity (top50)\nPearson r={pearson_r:.3f}, Spearman r={spearman_r:.3f}')
plt.grid(True, alpha=0.3)

# 标注几个极端点
for _, row in df.nlargest(3, 'rank_diff').iterrows():
    plt.annotate(row['chembl_id'], (row['docking_score'], row['predicted_pIC50']),
                 fontsize=7, alpha=0.8)

plt.tight_layout()
plt.savefig('correlation_plot.png', dpi=150)
print("\n散点图已保存到 correlation_plot.png")

# 保存带排名差异的表
df.to_csv('top50_with_ranks.csv', index=False)
print("完整结果已保存到 top50_with_ranks.csv")

python correlation_analysis.py

# 相关系数：
# # Pearson r > 0.3：对接打分和 QSAR 预测有中等正相关，说明两者大体一致
# # Pearson r 接近 0：两者基本无关，说明对接打分函数和 QSAR 模型关注的结构特征不同
# # Pearson r < 0：两者矛盾，对接打分好的 QSAR 预测反而低，需要警惕

### 简历：对 top50 候选分子进行对接打分与 QSAR 预测活性的相关性分析，发现两者在截断样本上无显著正相关（Spearman r = -0.225），进一步识别出对接打分与 QSAR 预测分歧最大的分子（如 CHEMBL5571836），作为后续 MD 验证的优先对象。
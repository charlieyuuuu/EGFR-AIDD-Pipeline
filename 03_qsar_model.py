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

print("1. 加载数据...")
df = pd.read_csv('CHEMBL203_actives.csv')
print(f"   原始数据: {len(df)} 条")

df = df.dropna(subset=['smiles', 'pChEMBL'])
df['pChEMBL'] = pd.to_numeric(df['pChEMBL'], errors='coerce')
df = df.dropna(subset=['pChEMBL'])
print(f"   清理后: {len(df)} 条")


print("2. 生成 Morgan 指纹 (ECFP4)...")
def smiles_to_fingerprint(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return np.array(fp)


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


print("3. 划分训练集/测试集并训练随机森林...")
# 80% 训练，20% 测试
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)
print(f"   训练集: {len(X_train)} 条, 测试集: {len(X_test)} 条")


model = RandomForestRegressor(
    n_estimators=100,
    max_depth=20,
    min_samples_leaf=3,
    random_state=42,
    n_jobs=-1
)
model.fit(X_train, y_train)


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


print("5. 预测 top50 对接结果的活性...")

try:
    top50 = pd.read_csv('top50_docking.csv')
except FileNotFoundError:
    top50 = pd.read_csv('docking_results.csv').sort_values('affinity').head(50)
    top50.to_csv('top50_docking.csv', index=False)


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


df = pd.read_csv('top50_qsar_prediction.csv')
print(f"化合物数: {len(df)}")

df['docking_score'] = -df['affinity']


pearson_r, pearson_p = pearsonr(df['docking_score'], df['predicted_pIC50'])
spearman_r, spearman_p = spearmanr(df['docking_score'], df['predicted_pIC50'])

print(f"\nPearson 相关系数: {pearson_r:.3f} (p={pearson_p:.4f})")
print(f"Spearman 相关系数: {spearman_r:.3f} (p={spearman_p:.4f})")
print("\n解释: 相关系数 > 0 表示对接打分越好，QSAR 预测活性也越高")


df['rank_docking'] = df['docking_score'].rank(ascending=False)
df['rank_qsar'] = df['predicted_pIC50'].rank(ascending=False)
df['rank_diff'] = abs(df['rank_docking'] - df['rank_qsar'])

print("\n排名差异最大的 10 个（对接和 QSAR 分歧最大）:")
print(df.nlargest(10, 'rank_diff')[
    ['chembl_id', 'affinity', 'predicted_pIC50', 'rank_docking', 'rank_qsar']
].to_string(index=False))

plt.figure(figsize=(8, 6))
plt.scatter(df['docking_score'], df['predicted_pIC50'], alpha=0.7, edgecolors='k')
plt.xlabel('Docking Score (-affinity, kcal/mol)')
plt.ylabel('Predicted pIC50')
plt.title(f'Docking Score vs QSAR Predicted Activity (top50)\nPearson r={pearson_r:.3f}, Spearman r={spearman_r:.3f}')
plt.grid(True, alpha=0.3)

for _, row in df.nlargest(3, 'rank_diff').iterrows():
    plt.annotate(row['chembl_id'], (row['docking_score'], row['predicted_pIC50']),
                 fontsize=7, alpha=0.8)

plt.tight_layout()
plt.savefig('correlation_plot.png', dpi=150)
print("\n散点图已保存到 correlation_plot.png")

df.to_csv('top50_with_ranks.csv', index=False)
print("完整结果已保存到 top50_with_ranks.csv")

python correlation_analysis.py

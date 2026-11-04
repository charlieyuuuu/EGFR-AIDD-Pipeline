"""
01_data_prep.py
Description: Filter ChEMBL compounds using Lipinski's rules and prepare 3D PDBQT ligands.
"""
# 从ChEMBL下载小分子数据库
# 小分子蛋白是EGFR
# 写一个下载脚本，nano download_chembl.py，内容可以根据你的需求调整。下面这个例子是按靶点下载所有有活性数据的化合物，非常适合做虚拟筛选的练习集：
from chembl_webresource_client.new_client import new_client
import pandas as pd
from rdkit import Chem

#  UniProt ID
target_uniprot = "P00533"

# 1. 找到靶点对应的 ChEMBL 内部 ID
target = new_client.target
target_query = target.search(target_uniprot).filter(
    organism__iexact='Homo sapiens'
).only(['target_chembl_id', 'pref_name'])
target_info = target_query[0]
target_chembl_id = target_info['target_chembl_id']
print(f"靶点: {target_info['pref_name']} ({target_chembl_id})")

# 2. 拉取该靶点的所有活性数据
activity = new_client.activity
res = activity.filter(
    target_chembl_id=target_chembl_id,
    standard_type="IC50",  # 也可以换成 Kd, Ki 等
    standard_relation="=",
    pchembl_value__isnull=False
).only([
    'molecule_chembl_id',
    'canonical_smiles',
    'standard_value',
    'standard_units',
    'pchembl_value'
])

# 3. 去重并整理成 DataFrame
records = []
seen = set()
for r in res:
    cid = r['molecule_chembl_id']
    if cid in seen:
        continue
    seen.add(cid)
    records.append({
        'chembl_id': cid,
        'smiles': r['canonical_smiles'],
        'IC50': r['standard_value'],
        'units': r['standard_units'],
        'pChEMBL': r['pchembl_value']
    })

df = pd.DataFrame(records)
print(f"共获取 {len(df)} 个唯一化合物")
df.to_csv(f'{target_uniprot}_actives.csv', index=False)
print(f"已保存到 {target_uniprot}_actives.csv")

python download_chembl.py

# 用 RDKit 做预处理和 Lipinski 过滤
nano filter_chembl.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

df = pd.read_csv('CHEMBL203_actives.csv')
print(f"原始数据: {len(df)} 条")

records = []
invalid = 0

for _, row in df.iterrows():
    smi = row['smiles']
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        invalid += 1
        continue

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1

    records.append({
        'chembl_id': row['chembl_id'],
        'smiles': smi,
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'HBD': hbd,
        'HBA': hba,
        'violations': violations,
        'pChEMBL': row['pChEMBL'],
    })

result = pd.DataFrame(records)
print(f"有效分子: {len(result)} 条")
print(f"无效 SMILES: {invalid} 条")

# 保留违规数 <= 1 的分子
filtered = result[result['violations'] <= 1].copy()
print(f"通过 Lipinski 过滤 (violations <= 1): {len(filtered)} 条")

filtered.to_csv('CHEMBL203_druglike.csv', index=False)
print("已保存到 CHEMBL203_druglike.csv")

# 看一眼分布
print("\nMW 范围:", filtered['MW'].min(), "-", filtered['MW'].max())
print("LogP 范围:", filtered['LogP'].min(), "-", filtered['LogP'].max())

python filter_chembl.py

# 构建抽样 + 批量生成 3D 构象
nano sample_and_embed.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os

df = pd.read_csv('CHEMBL203_druglike.csv')

# 排除 MW > 600 的极端分子
df = df[df['MW'] <= 600].copy()
print(f"排除 MW > 600 后: {len(df)} 条")

# 随机抽 100 个，固定随机种子保证可复现
sample = df.sample(n=100, random_state=42).reset_index(drop=True)

os.makedirs('ligands_3d', exist_ok=True)
os.makedirs('ligands_pdbqt', exist_ok=True)

success = 0
failed = []

for i, row in sample.iterrows():
    cid = row['chembl_id']
    smi = row['smiles']
    mol = Chem.MolFromSmiles(smi)

    if mol is None:
        failed.append((cid, 'SMILES parse failed'))
        continue

    mol = Chem.AddHs(mol)

    # 生成 3D 构象
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)

    if status != 0:
        failed.append((cid, '3D embedding failed'))
        continue

    # 力场优化
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        failed.append((cid, f'MMFF failed: {e}'))
        continue

    # 保存 SDF
    sdf_path = f'ligands_3d/{cid}.sdf'
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    success += 1
    if (i + 1) % 20 == 0:
        print(f"已处理 {i+1}/100，成功 {success}")

print(f"\n完成: 成功生成 3D 构象 {success}/100")
if failed:
    print(f"失败 {len(failed)} 个:")
    for cid, reason in failed[:10]:
        print(f"  {cid}: {reason}")

sample.to_csv('sample_100.csv', index=False)
print("\n抽样列表已保存到 sample_100.csv")

python sample_and_embed.py

#排除 MW > 600：避免大分子在 3D 生成时失败
#ETKDGv3 + MMFFOptimizeMolecule：这是 RDKit 生成 3D 构象的标准组合，ETKDG 负责生成初始坐标，MMFF 负责力场优化
#逐个保存 SDF：每个分子单独一个文件，方便后续转 PDBQT 和批量对接
#记录失败原因：有些分子可能因为环张力或特殊结构无法生成 3D，记录下来，不影响整体流程


# 将所有的小分子转成3D构象
nano embed_all.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import time

df = pd.read_csv('CHEMBL203_druglike.csv')

# 排除 MW > 600 的极端分子
df = df[df['MW'] <= 600].copy()
print(f"待处理分子总数: {len(df)}")

os.makedirs('ligands_3d', exist_ok=True)

success = 0
skipped = 0
failed = []

start_time = time.time()

for i, row in df.iterrows():
    cid = row['chembl_id']
    smi = row['smiles']
    sdf_path = f'ligands_3d/{cid}.sdf'

    # 已经生成过就跳过
    if os.path.exists(sdf_path):
        skipped += 1
        continue

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        failed.append((cid, 'SMILES parse failed'))
        continue

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)

    if status != 0:
        failed.append((cid, '3D embedding failed'))
        continue

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        failed.append((cid, f'MMFF failed: {e}'))
        continue

    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    success += 1

    if (success + skipped) % 500 == 0:
        elapsed = time.time() - start_time
        print(f"进度: {success + skipped}/{len(df)}  "
              f"成功 {success}  跳过 {skipped}  失败 {len(failed)}  "
              f"耗时 {elapsed:.0f}s")

print(f"\n完成")
print(f"新生成: {success}")
print(f"已跳过: {skipped}")
print(f"失败: {len(failed)}")
print(f"总耗时: {time.time() - start_time:.0f}s")

if failed:
    pd.DataFrame(failed, columns=['chembl_id', 'reason']).to_csv(
        'embed_failed.csv', index=False
    )
    print("失败列表已保存到 embed_failed.csv")

nohup python embed_all.py > embed_all.log 2>&1 &
# nohup：忽略挂断信号，SSH 断开也不影响
# # > embed_all.log 2>&1：把标准输出和错误都写到日志文件
# # &：放到后台运行
# # 运行后会返回一个进程号，比如 [1] 133800。

#查看进度
tail -f embed_all.log

#看完按 Ctrl+C 退出 tail，脚本本身不受影响，继续在后台跑。
#或者不进入实时模式，直接看最后几行：
tail -20 embed_all.log


# 批量转 PDBQT
# Meeko 是 AutoDock 官方推荐的配体准备工具，专门为 Vina 设计，对原子类型、电荷和扭转树（torsion tree）的处理更可靠
python -m pip install meeko

# Meeko 依赖的 scipy 没装。这是 Meeko 的 receptor_pdbqt 模块需要 scipy，但 pip 安装时不会自动拉取。
python -m pip install scipy

# gemmi是 Meeko 的另一个依赖，专门用来处理大分子结构。
python -m pip install gemmi

# 装完验证一下：
python -c "import meeko; print('Meeko OK')"

nano convert_to_pdbqt.py

import os
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

input_dir = 'ligands_3d'
output_dir = 'ligands_pdbqt'
os.makedirs(output_dir, exist_ok=True)

sdf_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.sdf')])
print(f"待转换 SDF 文件数: {len(sdf_files)}")

success = 0
skipped = 0
failed = []

for i, fname in enumerate(sdf_files):
    cid = fname.replace('.sdf', '')
    out_path = os.path.join(output_dir, f'{cid}.pdbqt')

    # 跳过已转换的
    if os.path.exists(out_path):
        skipped += 1
        continue

    sdf_path = os.path.join(input_dir, fname)
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = suppl[0] if len(suppl) > 0 else None

    if mol is None:
        failed.append((cid, 'RDKit read failed'))
        continue

    try:
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])

        if not is_ok:
            failed.append((cid, f'Meeko: {err}'))
            continue

        with open(out_path, 'w') as f:
            f.write(pdbqt_string)

        success += 1

    except Exception as e:
        failed.append((cid, str(e)))

    if (i + 1) % 500 == 0:
        print(f"进度: {i+1}/{len(sdf_files)}  成功 {success}  跳过 {skipped}  失败 {len(failed)}")

print(f"\n完成")
print(f"新转换: {success}")
print(f"已跳过: {skipped}")
print(f"失败: {len(failed)}")

if failed:
    import pandas as pd
    pd.DataFrame(failed, columns=['chembl_id', 'reason']).to_csv(
        'pdbqt_failed.csv', index=False
    )
    print("失败列表已保存到 pdbqt_failed.csv")

nohup python convert_to_pdbqt.py > convert_pdbqt.log 2>&1 &
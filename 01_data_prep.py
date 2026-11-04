"""
01_data_prep.py
Description: Filter ChEMBL compounds using Lipinski's rules and prepare 3D PDBQT ligands.
"""
# 从ChEMBL下载小分子数据库
# 小分子蛋白是EGFR
from chembl_webresource_client.new_client import new_client
import pandas as pd
from rdkit import Chem

#  UniProt ID
target_uniprot = "P00533"

target = new_client.target
target_query = target.search(target_uniprot).filter(
    organism__iexact='Homo sapiens'
).only(['target_chembl_id', 'pref_name'])
target_info = target_query[0]
target_chembl_id = target_info['target_chembl_id']
print(f"靶点: {target_info['pref_name']} ({target_chembl_id})")

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

filtered = result[result['violations'] <= 1].copy()
print(f"通过 Lipinski 过滤 (violations <= 1): {len(filtered)} 条")

filtered.to_csv('CHEMBL203_druglike.csv', index=False)
print("已保存到 CHEMBL203_druglike.csv")

print("\nMW 范围:", filtered['MW'].min(), "-", filtered['MW'].max())
print("LogP 范围:", filtered['LogP'].min(), "-", filtered['LogP'].max())

python filter_chembl.py

nano sample_and_embed.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os

df = pd.read_csv('CHEMBL203_druglike.csv')

df = df[df['MW'] <= 600].copy()
print(f"排除 MW > 600 后: {len(df)} 条")

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



nano embed_all.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import time

df = pd.read_csv('CHEMBL203_druglike.csv')

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
tail -f embed_all.log
tail -20 embed_all.log


python -m pip install meeko
python -m pip install scipy
python -m pip install gemmi
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

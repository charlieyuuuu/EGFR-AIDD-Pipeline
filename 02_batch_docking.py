"""
02_batch_docking.py
Description: Automated batch docking using AutoDock Vina via subprocess.
"""
# 直接用 Meeko 的命令行工具将pdb格式转换成pdbqt：
# 先安装一下Meeko依赖包prody
python -m pip install prody
python -c "import prody; print('prody OK')"

mk_prepare_receptor.py -i EGFR.pdb -o EGFR_receptor -p
# -i：输入的 PDB 文件
# # -o：输出文件前缀，会生成 EGFR_receptor.pdbqt
# # -p：保留配体坐标（通常受体准备需要去掉配体，但先跑通流程看它报不报错）

#第一步：确认受体 PDBQT 是“干净”的
#先快速看一眼文件，确认没有多余的分子被带进来：
grep "^HETATM" EGFR_receptor.pdbqt | head
grep "^ATOM" EGFR_receptor.pdbqt | wc -l

# 第二步：确定对接盒子（Grid Box）的位置
# 从原始 PDB 里提取共晶配体坐标
grep "^HETATM" 9BY4.pdb | awk '{print $4}' | sort -u
# Q6K 就是 9BY4 里与 EGFR 共晶的抑制剂

# 算 Q6K 的几何中心
nano get_ligand_center.py

import numpy as np

pdb_file = '9BY4.pdb'
ligand_resname = 'Q6K'

coords = []

with open(pdb_file, 'r') as f:
    for line in f:
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            if resname == ligand_resname:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

if coords:
    coords = np.array(coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)

    print(f"配体 {ligand_resname} 原子数: {len(coords)}")
    print(f"\n盒子中心坐标:")
    print(f"  --center_x {center[0]:.2f}")
    print(f"  --center_y {center[1]:.2f}")
    print(f"  --center_z {center[2]:.2f}")

    print(f"\n配体自身尺寸 (x, y, z): {span[0]:.2f}, {span[1]:.2f}, {span[2]:.2f}")
    print(f"建议盒子大小: 至少 {max(span)+10:.0f} Å")
else:
    print(f"没有找到残基 {ligand_resname}")


python get_ligand_center.py


# vina批量对接
nano batch_docking.py

import os
os.environ['PATH'] = '/home/zxsys/.local/bin:' + os.environ['PATH']
import subprocess
import pandas as pd
import time

VINA = '/home/zxsys/.local/bin/vina'

receptor = 'EGFR_receptor.pdbqt'
ligand_dir = 'ligands_pdbqt'
output_dir = 'docking_results'
os.makedirs(output_dir, exist_ok=True)

# 盒子参数
center_x, center_y, center_z = 50.62, -0.20, 18.68
size_x = size_y = size_z = 24
exhaustiveness = 8

ligands = sorted([f for f in os.listdir(ligand_dir) if f.endswith('.pdbqt')])
print(f"待对接配体总数: {len(ligands)}")

results = []
failed = []
start_time = time.time()

for i, lig in enumerate(ligands):
    cid = lig.replace('.pdbqt', '')
    out_path = os.path.join(output_dir, f'{cid}_out.pdbqt')

    # 跳过已完成的
    if os.path.exists(out_path):
        continue

    cmd = [
        VINA,
        '--receptor', receptor,
        '--ligand', os.path.join(ligand_dir, lig),
        '--center_x', str(center_x),
        '--center_y', str(center_y),
        '--center_z', str(center_z),
        '--size_x', str(size_x),
        '--size_y', str(size_y),
        '--size_z', str(size_z),
        '--exhaustiveness', str(exhaustiveness),
        '--out', out_path,
        '--cpu', '4',
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if proc.returncode != 0:
            failed.append((cid, f'vina return {proc.returncode}: {proc.stderr.strip()[:100]}'))
            continue

        # 从 stdout 解析最佳打分（Vina 1.2.5 不再写 --log 文件）
        affinity = None
        for line in proc.stdout.splitlines():
            line = line.strip()
            if line.startswith('1 '):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinity = float(parts[1])
                    except ValueError:
                        pass
                break

        if affinity is not None:
            results.append({'chembl_id': cid, 'affinity': affinity})
        else:
            failed.append((cid, 'no affinity parsed'))

    except subprocess.TimeoutExpired:
        failed.append((cid, 'timeout'))
    except Exception as e:
        failed.append((cid, str(e)))

    if (i + 1) % 100 == 0:
        elapsed = time.time() - start_time
        print(f"进度: {i+1}/{len(ligands)}  成功 {len(results)}  失败 {len(failed)}  耗时 {elapsed:.0f}s")

print(f"\n完成")
print(f"成功: {len(results)}")
print(f"失败: {len(failed)}")

if results:
    df = pd.DataFrame(results).sort_values('affinity')
    df.to_csv('docking_results.csv', index=False)
    print(f"结果已保存到 docking_results.csv")
    print(f"最佳打分: {df['affinity'].min():.2f} kcal/mol")

if failed:
    pd.DataFrame(failed, columns=['chembl_id', 'reason']).to_csv(
        'docking_failed.csv', index=False
    )
    print("失败列表已保存到 docking_failed.csv")


nohup python batch_docking.py > batch_docking.log 2>&1 &

tail -f batch_docking.log
ls docking_results/ | wc -l


## 存在小分子有多个构象的情况：
# “多个构象”指的是同一个分子在蛋白口袋里的不同结合姿势（binding poses）
# 要用vina_split进行分割，然后添加config文件进行打分：用 vina_split + --score_only 补上
nano rescore_missing.py

import os
import subprocess
import pandas as pd
import re

VINA = '/home/zxsys/.local/bin/vina'
VINA_SPLIT = '/home/zxsys/.local/bin/vina_split'
receptor = 'EGFR_receptor.pdbqt'

missing = pd.read_csv('missing_has_out.csv')['chembl_id'].tolist()
print(f"待补打分: {len(missing)}")

results = []
failed = []

for i, cid in enumerate(missing):
    out_path = f'docking_results/{cid}_out.pdbqt'
    if not os.path.exists(out_path):
        failed.append((cid, 'no output file'))
        continue

    # 拆分，取 pose 1
    try:
        subprocess.run(
            [VINA_SPLIT, '--input', out_path],
            capture_output=True, text=True, timeout=60
        )
    except Exception as e:
        failed.append((cid, f'split error: {e}'))
        continue

    pose1_path = f'docking_results/{cid}_out_ligand_1.pdbqt'
    if not os.path.exists(pose1_path):
        failed.append((cid, 'pose1 not generated'))
        continue

    # 对 pose1 打分（必须带盒子参数）
    score_cmd = [
        VINA, '--receptor', receptor,
        '--ligand', pose1_path,
        '--score_only',
        '--center_x', '50.62', '--center_y', '-0.20', '--center_z', '18.68',
        '--size_x', '24', '--size_y', '24', '--size_z', '24',
    ]
    try:
        score_proc = subprocess.run(
            score_cmd, capture_output=True, text=True, timeout=120
        )
    except Exception as e:
        failed.append((cid, f'score error: {e}'))
        os.remove(pose1_path)
        continue

    # 解析 "Estimated Free Energy of Binding : -8.485 (kcal/mol)"
    affinity = None
    for line in score_proc.stdout.splitlines():
        m = re.search(r'Estimated Free Energy of Binding\s*:\s*(-?\d+\.?\d*)', line)
        if m:
            affinity = float(m.group(1))
            break

    # 清理拆分文件
    if os.path.exists(pose1_path):
        os.remove(pose1_path)

    if affinity is not None:
        results.append({'chembl_id': cid, 'affinity': affinity})
    else:
        failed.append((cid, 'no affinity parsed'))

    if (i + 1) % 200 == 0:
        print(f"进度: {i+1}/{len(missing)}  成功 {len(results)}  失败 {len(failed)}")

print(f"\n完成: 成功 {len(results)}  失败 {len(failed)}")

if results:
    pd.DataFrame(results).to_csv('rescore_results.csv', index=False)
    print("已保存到 rescore_results.csv")

if failed:
    pd.DataFrame(failed, columns=['chembl_id', 'reason']).to_csv(
        'rescore_failed.csv', index=False
    )
    print("失败列表已保存到 rescore_failed.csv")

nohup python rescore_missing.py > rescore.log 2>&1 &

# 整合对接结果
nano merge.csv

import pandas as pd

df1 = pd.read_csv('docking_results.csv')
df2 = pd.read_csv('rescore_results.csv')
merged = pd.concat([df1, df2]).drop_duplicates(subset='chembl_id', keep='first')
merged = merged.sort_values('affinity')
merged.to_csv('docking_results_full.csv', index=False)

print(f"原始: {len(df1)}")
print(f"补分: {len(df2)}")
print(f"合并后: {len(merged)}")
print(f"最佳打分: {merged['affinity'].min():.2f}")
print(f"打分范围: {merged['affinity'].min():.2f} 到 {merged['affinity'].max():.2f}")

python merge.csv


# 对接结果排序，选 top50
nano top50.py

import pandas as pd

df = pd.read_csv('docking_results_full.csv')
top50 = df.sort_values('affinity').head(50)
top50.to_csv('top50_docking.csv', index=False)

print(f"总化合物: {len(df)}")
print(f"top50 打分范围: {top50['affinity'].min():.2f} 到 {top50['affinity'].max():.2f}")
print(top50.head(10).to_string(index=False))

python top50.py

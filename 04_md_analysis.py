"""
04_md_analysis.py
Description: Analyze GROMACS trajectories (.xvg) and plot SCI-standard RMSD curves.
"""
# 类型	chembl_id	对接打分	QSAR 预测 pIC50	选它的理由
# 共识分子	CHEMBL6001037	-10.65	8.92	对接和 QSAR 都排前列
# 共识分子	CHEMBL5920247	-10.55	8.89	同上，结构不同，增加多样性
# 对接极好但 QSAR 不认可	CHEMBL5571836	-11.57	6.01	全场对接最佳，QSAR 排名 40
# 对接极好但 QSAR 不认可	CHEMBL3742060	-10.72	5.55	对接排名 6.5，QSAR 排名 49，分歧最大之一
# QSAR 极好但对接一般	CHEMBL5959712	-10.11	8.62	对接排名 49，QSAR 排名 4
# QSAR 极好但对接一般	CHEMBL5805459	-10.23	8.63	对接排名 40，QSAR 排名 3


mk_export.py docking_results/CHEMBL6001037_out.pdbqt -s CHEMBL6001037_ligand.sdf
ls -la CHEMBL6001037_ligand.sdf
grep -c '^\$\$\$\$' CHEMBL6001037_ligand.sdf
nano extract_all_pose1.py

from rdkit import Chem
import os

candidates = [
    'CHEMBL6001037', 'CHEMBL5920247', 'CHEMBL5571836',
    'CHEMBL3742060', 'CHEMBL5959712', 'CHEMBL5805459'
]

for cid in candidates:
    sdf_path = f'{cid}_ligand.sdf'
    pdb_path = f'{cid}_pose1.pdb'

    if not os.path.exists(sdf_path):
        print(f"{cid}: SDF 不存在，跳过")
        continue

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = suppl[0] if len(suppl) > 0 else None
    if mol is None:
        print(f"{cid}: 读取失败")
        continue

    Chem.MolToPDBFile(mol, pdb_path)
    n_atoms = mol.GetNumAtoms()
    n_h = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'H')
    print(f"{cid}: {n_atoms} 原子, {n_h} 氢 -> {pdb_path}")

python extract_all_pose1.py


for cid in CHEMBL6001037 CHEMBL5920247 CHEMBL5571836 CHEMBL3742060 CHEMBL5959712 CHEMBL5805459; do
    echo -n "$cid 重复原子名: "
    awk '/^ATOM|^HETATM/ {print substr($0, 13, 4)}' ${cid}_pose1.pdb | sort | uniq -d | wc -l
done


conda activate acpype_new

nohup bash -c '
for cid in CHEMBL5920247 CHEMBL5571836 CHEMBL3742060 CHEMBL5959712 CHEMBL5805459; do
    echo "========== acpype: $cid =========="
    acpype -i ${cid}_pose1.pdb -b ${cid} -a gaff2 -n 0
done
' > acpype_batch.log 2>&1 &

for cid in CHEMBL6001037 CHEMBL5920247 CHEMBL5571836 CHEMBL3742060 CHEMBL5959712 CHEMBL5805459; do
    echo -n "$cid: "
    ls ${cid}.acpype/${cid}_GMX.itp 2>/dev/null && echo OK || echo MISSING
done


gmx pdb2gmx -f EGFR.pdb -o protein.gro -p topol.top -water tip3p -ff amber99sb-ildn -ignh

head -40 CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp
grep -n "moleculetype" CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp
cp CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp.bak
head -n 24 CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp > CHEMBL6001037_atomtypes.txt
tail -n +25 CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp > CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp.new
mv CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp.new CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp


nano topol.top


sed -i '27d' topol.top
sed -i '/#include "amber99sb-ildn.ff\/forcefield.itp"/a #include "CHEMBL6001037_atomtypes.txt"' topol.top
sed -i '/^\[ system \]/i #include "CHEMBL6001037.acpype\/CHEMBL6001037_GMX.itp"' topol.top
sed -i '48489a CHEMBL6001037       1' topol.top

grep -n "atomtypes\|CHEMBL6001037\|\[ system \]\|\[ molecules \]\|Protein_chain_A" topol.top


nano em.mdp

; em.mdp - energy minimization
integrator  = steep
emtol       = 100.0
emstep      = 0.01
nsteps      = 50000

; Neighbor searching
cutoff-scheme = Verlet
nstlist       = 10
ns_type       = grid
pbc           = xyz

; Electrostatics
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0

nano merge_gro.py

def merge_gro(protein_gro, ligand_gro, output_gro):
    with open(protein_gro) as f:
        plines = f.readlines()
    with open(ligand_gro) as f:
        llines = f.readlines()

    p_atoms = int(plines[1].strip())
    l_atoms = int(llines[1].strip())
    total = p_atoms + l_atoms

    out = [plines[0], f"{total}\n"] + plines[2:-1] + llines[2:-1] + [plines[-1]]

    with open(output_gro, 'w') as f:
        f.writelines(out)
    print(f"合并: {p_atoms} + {l_atoms} = {total} 原子 -> {output_gro}")

merge_gro('protein.gro', 'CHEMBL6001037.acpype/CHEMBL6001037_GMX.gro', 'complex_CHEMBL6001037.gro')

python merge_gro.py


## 进行能量最小化
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em2.tpr -maxwarn 2
gmx mdrun -v -deffnm em2 -ntomp 12

## 建盒子、加溶剂、加离子
# 定义盒子
gmx editconf -f em.gro -o box.gro -c -d 1.0 -bt dodecahedron

# 加溶剂
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# 加离子
nano ions.mdp

integrator  = steep
emtol       = 1000.0
nsteps      = 50000

cutoff-scheme = Verlet
nstlist       = 10
pbc           = xyz

coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0

gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
grep -A15 "\[ molecules \]" topol.top


gmx mdrun -v -deffnm em_pose2 -nb gpu -ntmpi 1 -ntomp 12
# NVT
gmx grompp -f nvt.mdp -c em_pose2.gro -p topol_pose2.top -o nvt_pose2.tpr -maxwarn 2
gmx mdrun -v -deffnm nvt_pose2 -nb gpu -ntmpi 1 -ntomp 8

# NPT
gmx grompp -f npt.mdp -c nvt_pose2.gro -p topol_pose2.top -o npt_pose2.tpr -maxwarn 2
gmx mdrun -v -deffnm npt_pose2 -nb gpu -ntmpi 1 -ntomp 8



### NVT 平衡
nano nvt.mdp

; nvt.mdp - NVT equilibration
integrator  = md
dt          = 0.002
nsteps      = 50000    ; 100 ps

nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

continuation = no
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4

cutoff-scheme = Verlet
nstlist     = 10
pbc         = xyz

coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300

pcoupl      = no

gen_vel     = yes
gen_temp    = 300
gen_seed    = -1

gmx grompp -f nvt.mdp -c solv_ions.gro -p topol.top -o nvt.tpr -maxwarn 2

gmx energy -f nvt.edr -o nvt_temp.xvg
grep -v "^[#@]" nvt_temp.xvg | awk '{sum+=$2; n++} END {print "平均温度:", sum/n, "K"}'



# NPT 平衡
nano npt.mdp

; npt.mdp - NPT equilibration
integrator  = md
dt          = 0.002
nsteps      = 50000    ; 100 ps

nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4

cutoff-scheme = Verlet
nstlist     = 10
pbc         = xyz

coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300

pcoupl      = C-rescale
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

gen_vel     = no

gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -deffnm npt -ntomp 12
echo "Density" | gmx energy -f npt.edr -o npt_density.xvg
echo "Pressure" | gmx energy -f npt.edr -o npt_pressure.xvg



nano md.mdp

; md.mdp - production run 10 ns
integrator  = md
dt          = 0.002
nsteps      = 5000000    ; 10 ns

nstxout-compressed = 5000   ; 每 10 ps 存一帧（10 ns 共 1000 帧）
nstenergy   = 5000
nstlog      = 5000

continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 1
lincs_order = 4

cutoff-scheme = Verlet
nstlist     = 10
pbc         = xyz

coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300

pcoupl      = C-rescale
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

gen_vel     = no

gmx grompp -f md.mdp -c npt2.gro -p topol.top -o md.tpr -maxwarn 2
nohup gmx mdrun -v -deffnm md -nb gpu -ntomp 8 > md.out 2>&1 &
tail -20 md.log

ls -la md.gro md.xtc md.edr


## RMSD 分析
# 1. 蛋白骨架 RMSD
gmx rms -s md.tpr -f md.xtc -o rmsd_backbone.xvg -tu ns -n index.ndx

# 2. 配体 RMSD（关键）
gmx rms -s md.tpr -f md.xtc -o rmsd_ligand.xvg -tu ns -n index.ndx

grep -v "^[#@]" rmsd_backbone.xvg | awk '{sum+=$2; n++} END {print "Backbone 平均 RMSD:", sum/n, "nm"}'
grep -v "^[#@]" rmsd_ligand.xvg | awk '{sum+=$2; n++; if($1>5) late+=$2; if($1>5) ln++} END {print "配体全程平均:", sum/n, "nm"; print "5-10ns 平均:", late/ln, "nm"}'
gmx select -s md.gro -select 'resname UNL' -on ligand.ndx
echo "0" | gmx trjconv -s md_pose2.tpr -f md_pose2_center.xtc -o final_pose2.pdb -dump 100000 -n index_pose2.ndx


## 生成核心动力学折线图数据
# 1. 蛋白骨架 RMSD 和 回旋半径 (Rg)
# 蛋白骨架 RMSD
echo -e "4\n4" | gmx rms -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o rmsd_backbone.xvg -tu ns

# 回旋半径 (蛋白没散架)
echo "1" | gmx gyrate -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o gyrate.xvg

# 2. 蛋白各残基的柔性 (RMSF)
echo "1" | gmx rmsf -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o rmsf.xvg -res


# 3. 氢键数量（用替代方案避开报错）
# 算配体和蛋白的极性接触（距离 < 0.35 nm 就等同于氢键）
# 提取蛋白极性原子和配体极性原子
gmx make_ndx -f md_pose2.tpr -o index_final.ndx

# 算最小距离
echo -e "20\n21" | gmx mindist -s md_pose2.tpr -f md_pose2_center.xtc -n index_final.ndx -od hbond_distance.xvg -tu ns

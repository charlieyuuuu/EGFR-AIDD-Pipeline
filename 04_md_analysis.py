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

# 最直接的替代方案：用 Meeko 导出 SDF
# 确保在 rdkit-env 环境下（Meeko 装在这里）
conda activate rdkit-env
cd ~/YM/aidd_learning

# 用 Meeko 导出 SDF
mk_export.py docking_results/CHEMBL6001037_out.pdbqt -s CHEMBL6001037_ligand.sdf

# 看文件是否生成
ls -la CHEMBL6001037_ligand.sdf

# 看里面有几个构象（应该和原 PDBQT 的 pose 数一致）
grep -c '^\$\$\$\$' CHEMBL6001037_ligand.sdf


# 提取 pose 1 并转成 acpype 能用的格式
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

# 检查 PDB 里的原子名是否唯一
for cid in CHEMBL6001037 CHEMBL5920247 CHEMBL5571836 CHEMBL3742060 CHEMBL5959712 CHEMBL5805459; do
    echo -n "$cid 重复原子名: "
    awk '/^ATOM|^HETATM/ {print substr($0, 13, 4)}' ${cid}_pose1.pdb | sort | uniq -d | wc -l
done
# 如果没有输出，说明原子名唯一，可以直接进 acpype。如果有输出，说明 Meeko 导出的 PDB 里原子名有重复（这在某些配体上会发生），acpype 会报错。

# 将小分子转化成拓扑格式
cd ~/YM/aidd_learning
conda activate acpype_new

nohup bash -c '
for cid in CHEMBL5920247 CHEMBL5571836 CHEMBL3742060 CHEMBL5959712 CHEMBL5805459; do
    echo "========== acpype: $cid =========="
    acpype -i ${cid}_pose1.pdb -b ${cid} -a gaff2 -n 0
done
' > acpype_batch.log 2>&1 &

# 查 6 个配体的拓扑文件：
for cid in CHEMBL6001037 CHEMBL5920247 CHEMBL5571836 CHEMBL3742060 CHEMBL5959712 CHEMBL5805459; do
    echo -n "$cid: "
    ls ${cid}.acpype/${cid}_GMX.itp 2>/dev/null && echo OK || echo MISSING
done


### 开始将蛋白格式转化
# acpype 生成的是 GAFF2（General AMBER Force Field）拓扑。GROMACS 里必须用 AMBER 家族的蛋白力场（如 amber99sb-ildn 或 amber14sb）来配，不能用 CHARMM 或 OPLS。混用力场会导致原子类型冲突，grompp 直接失败。
gmx pdb2gmx -f EGFR.pdb -o protein.gro -p topol.top -water tip3p -ff amber99sb-ildn -ignh
# 参数说明：
# # -ff amber99sb-ildn：选 AMBER 力场，和配体的 GAFF2 兼容
# # -ignh：忽略 PDB 里的氢，让 pdb2gmx 按力场规则重新加氢（推荐，因为晶体结构的氢位置不一定符合力场）
# # -water tip3p：水模型，后续加溶剂要用一致的
# # 跑完后确认：protein.gro 和 topol.top 生成了。topol.top 里应该有 #include "amber99sb-ildn.ff/forcefield.itp"

# 处理配体 ITP（去掉头尾，只留 moleculetype）
# 看ITP结构：
head -40 CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp

# 接下来要做的是把 ITP 拆成两部分：[ atomtypes ] 合并到 topol.top，[ moleculetype ] 及之后的内容保留在 ITP 里。然后修改 topol.top 的 include 和 [ molecules ]。
# 第一步：确认 [ moleculetype ] 在第几行
grep -n "moleculetype" CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp

# 拆分 ITP
cd ~/YM/aidd_learning
# 备份
cp CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp.bak
# 提取 atomtypes（第 1-24 行）
head -n 24 CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp > CHEMBL6001037_atomtypes.txt
# 保留 moleculetype 及之后（从第 25 行开始）
tail -n +25 CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp > CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp.new
mv CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp.new CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp

## 合并配体 ITP 到 topol.top
# 打开 topol.top：
nano topol.top

# 执行三处修改
# 删掉第27行
sed -i '27d' topol.top
# 第一处：在 forcefield include 之后加 atomtypes
sed -i '/#include "amber99sb-ildn.ff\/forcefield.itp"/a #include "CHEMBL6001037_atomtypes.txt"' topol.top
# 第二处：在 [ system ] 之前加配体 ITP include
sed -i '/^\[ system \]/i #include "CHEMBL6001037.acpype\/CHEMBL6001037_GMX.itp"' topol.top
# 第三处：在 [ molecules ] 里加配体
sed -i '48489a CHEMBL6001037       1' topol.top

# 🔍 修改后确认
grep -n "atomtypes\|CHEMBL6001037\|\[ system \]\|\[ molecules \]\|Protein_chain_A" topol.top
# 应该看到：
## #include "CHEMBL6001037_atomtypes.txt" 在 forcefield include 之后
## #include "CHEMBL6001037.acpype/CHEMBL6001037_GMX.itp" 在 [ system ] 之前
## [ molecules ] 里有 Protein_chain_A 1 和 CHEMBL6001037 1

### 拓扑合并完成
#现在 topol.top 的结构是：
# 1. forcefield include
# 2. 配体 atomtypes include
# 3. 蛋白 moleculetype + atoms/bonds/dihedrals 等
# 4. 配体 ITP include（在 [ system ] 之前）
# 5. [ molecules ] 里有蛋白和配体各 1 个


# 下一步：合并坐标（complex.gro）
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


# 用 Python 把 protein.gro 和 CHEMBL6001037.acpype/CHEMBL6001037_GMX.gro 合并：
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
# -nb gpu 让非键相互作用在 GPU 上算，加速明显

## 建盒子、加溶剂、加离子
# 定义盒子
gmx editconf -f em.gro -o box.gro -c -d 1.0 -bt dodecahedron
# -d 1.0 表示蛋白到盒子边缘至少 1.0 nm，-bt cubic 是立方盒子。对于球形蛋白，也可以选 dodecahedron（菱形十二面体），体积更小、更省水，但立方盒子更简单。
# 加溶剂
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
# 这会往盒子里填水分子，并自动在 topol.top 的 [ molecules ] 里加一行 SOL。注意：加水后 topol.top 会被修改，你之前的蛋白和配体行还在，只是多了 SOL。
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
# 这会生成 ions.tpr。然后用 genion 替换水分子为离子：
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
grep -A15 "\[ molecules \]" topol.top


## pose2参数
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
# pose1 gmx mdrun -v -deffnm nvt -nb gpu -ntmpi 1 -ntomp 12

# 先检查温度是否稳定
gmx energy -f nvt.edr -o nvt_temp.xvg
# 交互式提示时输入 Temperature 和 0（退出），然后看 nvt_temp.xvg 里的平均值：
grep -v "^[#@]" nvt_temp.xvg | awk '{sum+=$2; n++} END {print "平均温度:", sum/n, "K"}'
# 应该在 300 K 附近（±5 K 以内算正常）。如果偏差很大（比如 > 320 K 或 < 280 K），说明温度耦合有问题，需要检查 nvt.mdp 里的 tc-grps 和 tau_t。


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
# continuation = yes：从 NVT 继续，不重新生成速度
# # gen_vel = no：NPT 不需要重新初始化速度
# # pcoupl = C-rescale：GROMACS 2025 推荐的压力耦合方法，比 Berendsen 更准确
# # NPT 也是 100 ps，纯 CPU 大约 50 分钟。跑完后检查密度是否稳定在 1000 kg/m³ 左右（水）以及压力是否在 1 bar 附近波动。
# 密度（水应该在 1000 kg/m³ 左右）
echo "Density" | gmx energy -f npt.edr -o npt_density.xvg
# 压力（应该在 1 bar 附近波动）
echo "Pressure" | gmx energy -f npt.edr -o npt_pressure.xvg
# 判断标准
# 密度：如果体系里有蛋白和配体，整体密度会略高于纯水（约 1010-1050 kg/m³ 都是正常的）。如果密度在 1000 左右且波动小，说明 NPT 平衡成功。
# 压力：瞬时压力波动很大（几百 bar 都正常），但平均值应该在 1 bar 附近。如果平均值偏差很大（比如 > 100 bar），说明压力耦合可能有问题。


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
nvidia-smi

#  输出文件：
ls -la md.gro md.xtc md.edr
# md.xtc 是轨迹文件（约 1000 帧），md.gro 是最终结构，md.edr 是能量文件。


## 立刻做 RMSD 分析
# 1. 蛋白骨架 RMSD
gmx rms -s md.tpr -f md.xtc -o rmsd_backbone.xvg -tu ns -n index.ndx
# 交互式提示时，按顺序输入：
# # 两次都选择（拟合组）：输入 4（Backbone）

# 2. 配体 RMSD（关键）
gmx rms -s md.tpr -f md.xtc -o rmsd_ligand.xvg -tu ns -n index.ndx
# 交互式提示时，按顺序输入：
# # 第一次选择（拟合组）：输入 4（Backbone）
# # 第二次选择（计算组）：输入 13（UNL）
# 注意：UNL 是配体的残基名（acpype 生成的配体默认叫 UNL）。如果 gmx rms 提示找不到 UNL，可以用 gmx select 或 gmx make_ndx 确认配体组名

# 查看结果
# 蛋白骨架 RMSD 平均值
grep -v "^[#@]" rmsd_backbone.xvg | awk '{sum+=$2; n++} END {print "Backbone 平均 RMSD:", sum/n, "nm"}'
# 配体 RMSD 平均值和后期趋势
grep -v "^[#@]" rmsd_ligand.xvg | awk '{sum+=$2; n++; if($1>5) late+=$2; if($1>5) ln++} END {print "配体全程平均:", sum/n, "nm"; print "5-10ns 平均:", late/ln, "nm"}'
#注意 RMSD 单位是 nm（GROMACS 默认），1 nm = 10 Å。配体 RMSD 如果在 0.1-0.2 nm（1-2 Å）内稳定，说明结合良好

# # 看配体是否还在蛋白附近
gmx select -s md.gro -select 'resname UNL' -on ligand.ndx


# 核心目标是证明这个分子在 100 ns 里稳定结合。现在你手头的证据已经极其充分了：
# 1. 配体 RMSD 稳定在 0.3-0.6 nm（证明没解离）。
# 2. 最小距离稳定在 1.7 Å（证明紧密贴合，甚至有氢键）
# 漂亮的图，直接去提取第 90 ns 或 100 ns 的 PDB 结构，在 PyMOL 里画一下氢键，比 gmx hbond 直观一万倍：
# 提取最后 1 帧（100 ns）
echo "0" | gmx trjconv -s md_pose2.tpr -f md_pose2_center.xtc -o final_pose2.pdb -dump 100000 -n index_pose2.ndx
# 在 PyMOL 里打开 final_pose2.pdb，用 distance hbonds, resn UNL, protein, 3.5, mode=2，截图放进报告里。


## 第一弹：生成核心动力学折线图数据
# 你之前跑的 RMSD 是局部的，现在我们用 PBC 修正后的 md_pose2_center.xtc 来生成全局的曲线。
# 1. 蛋白骨架 RMSD 和 回旋半径 (Rg)
# 蛋白骨架 RMSD
echo -e "4\n4" | gmx rms -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o rmsd_backbone.xvg -tu ns

# 回旋半径 (证明蛋白没散架)
echo "1" | gmx gyrate -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o gyrate.xvg

# 2. 蛋白各残基的柔性 (RMSF)
# 算 RMSF 需要先做一下平均，不然图形会很乱：
echo "1" | gmx rmsf -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o rmsf.xvg -res


# 3. 氢键数量（用替代方案避开报错）
# 既然 gmx hbond 不支持 GAFF2，我们直接算配体和蛋白的极性接触（只要距离 < 0.35 nm 就等同于氢键）
# 提取蛋白极性原子和配体极性原子
gmx make_ndx -f md_pose2.tpr -o index_final.ndx

1 & a O* N*
name 20 Protein_Polar
13 & a O* N*
name 21 Ligand_Polar
q
# 用这两个组算最小距离
echo -e "20\n21" | gmx mindist -s md_pose2.tpr -f md_pose2_center.xtc -n index_final.ndx -od hbond_distance.xvg -tu ns
# 画图提示：拿到 .xvg 文件后，把数据导入 Origin，横坐标设为 0-100 ns，纵坐标用 nm 或 Å。你的目标图应该是：RMSD 平稳、Rg 平稳、RMSF 在结合口袋（残基 700-800）有一个明显的低谷。


## 第二弹：提取漂亮的 3D 结合模式图
# 1. 从轨迹中心提取代表结构（比如第 90 ns）
# 剥离水、离子，只保留蛋白和配体
echo -e "1 | 13" | gmx trjconv -s md_pose2.tpr -f md_pose2_center.xtc -n index_pose2.ndx -o complex_dry.pdb -dump 90000 -pbc mol -center

# 2. 在 PyMOL 里美化出图
# 用 PyMOL 打开 complex_dry.pdb，然后在 PyMOL 命令行（不是 Linux 终端）输入以下命令：[ray 的 PyMOL 渲染插件，或者直接用自带的 ray 命令，出来的图绝对有 SCI 范儿。]
# 1. 设置为白色背景
bg_color white
# 2. 隐藏氢原子 (让图变干净)
hide everything, hydrogens
# 3. 美化蛋白 (卡通模型，浅青色)
show cartoon, polymer
color palecyan, polymer
# 4. 美化配体 (绿色棍状)
show sticks, resn UNL
color green, resn UNL
# 5. 显示口袋残基 (距离配体 4 埃以内的残基)
select pocket, byres polymer within 4.0 of resn UNL
show sticks, pocket
color cyan, pocket
# 6. 自动寻找并画出氢键 (黄色虚线)
distance hbonds, resn UNL, polymer, 3.5, mode=2
show dashes, hbonds
color yellow, hbonds
# 7. 设置光线追踪，输出高清图
ray 1200, 900
png binding_mode_pose2.png
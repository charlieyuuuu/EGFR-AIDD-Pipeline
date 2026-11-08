# EGFR 靶点新型抑制剂的 AI 辅助全栈虚拟筛选与动力学验证
# AI-Assisted Full-Stack Virtual Screening & MD Validation Pipeline for EGFR Inhibitors

## 📌 项目概述 / Project Overview
本项目展示了一个端到端（End-to-End）的 AIDD/CADD 研发流水线，用于发现靶向 EGFR 的新型小分子抑制剂。工作流涵盖了**基于结构的虚拟筛选**、**机器学习（QSAR）交叉验证**、**100 纳秒（ns）全原子分子动力学（MD）模拟**、**MM-PBSA 结合自由能计算**、**ADMET 早期毒性评估**以及**分子优先级排序与优化决策**。核心目的是通过物理化学与 AI 的多维交叉验证，剔除传统分子对接的假阳性结果，全面评估候选分子的成药潜力，并为后续湿实验指明方向。
This repository presents an end-to-end AIDD/CADD pipeline for discovering novel EGFR inhibitors. The workflow integrates **structure-based virtual screening**, **machine learning (QSAR) cross-validation**, **100 ns all-atom molecular dynamics (MD) simulations**, **MM-PBSA binding free energy calculations**, **ADMET evaluation**, and **compound prioritization & optimization strategies**. The primary goal is to eliminate false positives, evaluate candidates comprehensively, and guide subsequent wet-lab experiments.

## 🛠️ 技术栈 / Tech Stack
*   **机器学习与化学信息学 (AIDD) / AI & Cheminformatics:** Scikit-learn (随机森林/Random Forest), RDKit, Pandas (Morgan 分子指纹/Fingerprints), ADMET-AI
*   **结构生物学与虚拟筛选 (CADD) / Structure-Based Design:** AutoDock Vina (v1.2.5), Meeko, OpenBabel
*   **分子动力学与热力学 (MD) / Molecular Dynamics:** GROMACS 2018 (GPU 加速/accelerated), acpype (GAFF2 力场/force field), gmx_MMPBSA
*   **自动化与可视化 / Automation & Visualization:** Linux Shell, Python `subprocess`, Matplotlib, PyMOL

## 📂 研发流水线概述 / Pipeline Overview
1.  **靶点结构选择与准备 / Target Structure Selection**：基于 PDB 数据库对 EGFR 激酶域结构进行系统评估。综合考虑**分辨率（9BY4 为 2.31 Å，满足 ≤ 2.5 Å 的优选标准）**、**共晶配体存在（含非共价抑制剂，明确界定 ATP 结合口袋边界）** 及 **野生型状态（无耐药突变）**，最终选用 **PDB: 9BY4** 用于分子对接，保证了格点设置的准确性与初始结构的可靠性。
    Systematically evaluated EGFR structures from the PDB. Considering **resolution (9BY4 at 2.31 Å)**, the presence of a **co-crystal ligand**, and **wild-type status**, **PDB: 9BY4** was selected to ensure accurate grid box setting.
2.  **数据准备 / Data Preparation**：从 ChEMBL 数据库中提取 10,734 个具有 EGFR 活性的化合物。利用 RDKit 进行 Lipinski 类药性五规则过滤，并生成 3D 构象及 PDBQT 格式文件。
    Filtered 10,734 EGFR-active compounds from the ChEMBL database using Lipinski's Rule of Five via RDKit. Generated 3D conformers and PDBQT formats.
3.  **高通量虚拟筛选与 Redocking 验证 (CADD) / High-Throughput Virtual Screening & Validation**：为确保对接参数的可靠性，本研究首先进行了 **Redocking 方法学验证**。通过将 9BY4 的共晶配体（Q6K）重新对接至原活性口袋，计算得出最优构象与晶体结构的 **RMSD 为 1.487 Å（< 2.0 Å）**，证实了对接参数（Exhaustiveness=16, 24Å 盒子）的精准性。随后，在 Linux 服务器上编写 Python 自动化脚本，调用 AutoDock Vina 对约 7,500 个化合物完成批量化对接，并根据结合能初筛出 Top 50 候选分子。
    Prior to screening, **Redocking validation** was performed to ensure reliability. By re-docking the co-crystal ligand (Q6K) of 9BY4, the **RMSD was calculated to be 1.487 Å (< 2.0 Å threshold)**, confirming the accuracy of docking parameters. Subsequently, automated the docking of ~7,500 compounds against the EGFR kinase domain using Python `subprocess` scripts, selecting the Top 50 candidates based on binding affinity.
4.  **AI / QSAR 交叉验证 (AIDD) / AI Cross-Validation**：提取 2048 位 ECFP4 Morgan 分子指纹，训练随机森林回归模型（测试集 R² = 0.682）预测 pIC50。利用该机器学习模型对 Top 50 对接结果进行活性交叉验证，识别假阳性。
    Extracted 2048-bit ECFP4 Morgan fingerprints and trained a Random Forest Regressor (Test R² = 0.682) to predict bioactivity. This ML model was used to cross-validate the Top 50 docking hits.
5.  **MD 体系搭建与模拟 (CADD) / MD Setup & Simulation**：针对 6 个代表性候选分子，利用 GROMACS 开展 100 ns 全原子动力学模拟。编写自动化脚本完成拓扑合并，并使用 GPU 大幅加速计算。
    Conducted 100 ns MD simulations on 6 representative complexes using GROMACS. Developed automated scripts for topology merging and utilized GPU acceleration.
6.  **MM-PBSA 自由能计算 / MM-PBSA Calculation**：基于 100 ns 轨迹提取 100 个代表构象，利用 `gmx_MMPBSA` 计算结合自由能并进行能量拆解。
    Calculated binding free energy using `gmx_MMPBSA` from 100 representative frames extracted from the trajectory, including energy decomposition.
7.  **ADMET 预测 / ADMET Prediction**：利用 ADMET-AI 工具包对候选分子进行 hERG 心脏毒性、Caco-2 渗透性及 CYP450 代谢预测，辅助早期筛选决策。
    Used ADMET-AI to predict hERG toxicity, Caco-2 permeability, and CYP450 metabolism for early-stage filtering.

## 🔬 核心科学发现 / Key Findings
*静态对接打分与真实的动力学行为存在显著差异，高分不等于可成药。通过结合 AI 模型、MD 模拟与热力学计算，我们揭示了以下核心现象：*
*Static docking scores and dynamic behavior differ significantly; a high score does not equal a good drug. By combining AI, MD, and thermodynamics, we revealed the following:*
*   **构象分化的假阳性陷阱 / Conformational Divergence**：Vina 打分较高（-10.65 kcal/mol）的 **CHEMBL6001037**，其主构象（pose 1）在 MD 模拟中表现出极不稳定的动力学特征（在 ~20 ns 首次脱离，~65 ns 完全水合解离）。而提取其次优构象（pose 2）重新模拟后，其通过口袋内的柔性适应（RMSD 0.4–0.65 nm），实现了全程 100 ns **未解离的稳定结合**。这证明了多构象采样与 MD 验证的必要性。
    **CHEMBL6001037** (pose 1) dissociated within ~20 ns. However, its secondary pose (pose 2) remained stably bound over 100 ns through flexible adaptation (RMSD 0.4–0.65 nm), proving the necessity of multi-pose sampling and MD validation.
*   **高分但非成药 / High Affinity but Poor Druggability**：**CHEMBL5571836** 虽然以 -11.57 kcal/mol 获得全场最高对接分，并在 MD（pose 1）中保持极高的结构刚性（~0.2 nm），但其 QSAR 活性预测不被认可，且 ADMET 预测显示极高的 hERG 毒性（0.90）与过快的 CYP3A4 代谢。因此，仅凭静态结合力无法判定其成药潜力。
    Despite the top score (-11.57 kcal/mol) and rigid MD stability, **CHEMBL5571836** was deprioritized due to poor QSAR prediction, high hERG toxicity (0.90), and fast metabolism.
*   **MM-PBSA 热力学验证 / Thermodynamic Validation**：对 Tier 1 候选分子（CHEMBL6001037 pose 2）进行结合自由能计算，得到 **ΔG ≈ -37 kcal/mol**。能量拆解表明，**范德华力（-49 kcal/mol）**与**静电力（-44 kcal/mol）**是复合物结合的主要驱动力。
    Binding free energy for the Tier 1 candidate was calculated as **ΔG ≈ -37 kcal/mol**, with **van der Waals (-49 kcal/mol)** and **electrostatic (-44 kcal/mol)** interactions as the main driving forces.

## 📊 结果展示 / Results

### 1. 方法学验证 (Redocking) / Methodological Validation
![Redocking Validation](results/Fig1_Redocking_Validation.png)
*图 1：Redocking 方法学验证。绿色为 9BY4 晶体原配体（Q6K），洋红色为 Vina 重对接最优构象，两者在原蛋白质结合口袋中高度重合（RMSD = 1.487 Å < 2.0 Å），证明了对接参数的极高可靠性。*
*Figure 1: Redocking validation. Green represents the crystal ligand (Q6K), and magenta represents the redocked pose, showing high overlap within the original binding pocket (RMSD = 1.487 Å < 2.0 Å).*

### 2. 候选复合物的动力学稳定性 / Dynamics Stability
![Ligand RMSD](results/Fig2_Ligand_RMSD.png)
*图 2：6 个候选复合物在 100 ns MD 模拟过程中的配体 RMSD 变化曲线。绿线（CHEMBL5571836 pose1）表现出极高的刚性（~0.2 nm）；蓝线（CHEMBL6001037 pose2）展现了口袋内 0.4–0.65 nm 的柔性适应，二者均未发生解离。*
*Figure 2: Ligand RMSD over 100 ns MD simulations. The green line shows high rigidity (~0.2 nm), while the blue line shows flexible adaptation (0.4–0.65 nm) without dissociation.*

### 3. 假阳性验证 / False Positive Validation
![Pose 1 Dissociation](results/Fig3_FalsePositive_Final.png)
*图 3：CHEMBL6001037 (pose 1) 的动力学解离过程。配体在 ~20 ns 首次脱离结合口袋，并在 ~65 ns 完全水合。曲线在 8–9 nm 与 0.5 nm 之间的剧烈跳跃为周期性边界（PBC）重入伪影，非重新结合。*
*Figure 3: Dissociation process of CHEMBL6001037 (pose 1). The ligand initially detached at ~20 ns and fully hydrated at ~65 ns. The rapid jumps between 8–9 nm and 0.5 nm are PBC re-entry artifacts, not re-binding.*

### 4. ADMET 早期成药性评估 / Early-stage ADMET Evaluation
我们利用 ADMET-AI 工具包对 Top 50 候选分子进行了评估。**关键结论**：结合能排名前列的候选分子虽然具备极佳的肠道吸收（Caco-2），但其 **hERG 心脏毒性预测值均偏高（> 0.8）**。这提示在后续先导化合物优化中，必须引入多参数优化（MPO）进行结构改造，以降低心脏毒性风险。
详细预测数据见：[`data/top50_ADMET_evaluated.csv`](data/top50_ADMET_evaluated.csv)

We utilized ADMET-AI to evaluate the Top 50 candidates. **Key conclusion**: Despite excellent intestinal absorption (Caco-2), the top-ranked molecules exhibit a **high risk of hERG cardiotoxicity (> 0.8)**. This highlights the necessity of Multi-Parameter Optimization (MPO) in subsequent lead optimization to mitigate cardiotoxicity risks.
Detailed prediction data: [`data/top50_ADMET_evaluated.csv`](data/top50_ADMET_evaluated.csv)

| ChEMBL ID | Affinity (kcal/mol) | hERG (Cardiotoxicity) | Caco2 (Absorption) | CYP3A4 (Metabolism) | Half-Life (hours) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **CHEMBL5571836** | -11.57 | 0.90 (High Risk) | -4.45 (Excellent) | 0.88 (Fast) | 15.62 |
| **CHEMBL337026** | -11.00 | 0.84 (High Risk) | -4.91 (Excellent) | 0.47 (Slow) | 347.91 |
| **CHEMBL338155** | -10.95 | 0.86 (High Risk) | -5.20 (Excellent) | 0.25 (Slow) | 377.84 |

## 💡 研发决策与结构优化建议 / Lead Optimization & Prioritization Strategy

综合分子对接打分、MD 动力学稳定性、MM-PBSA 热力学结合力与 ADMET 成药性，我们对候选化合物进行了多维度综合评估：
Based on virtual screening, MD stability, MM-PBSA thermodynamics, and ADMET predictions, we comprehensively evaluated the candidates:

### 1. 化合物综合优先级排序 / Compound Prioritization
*   **Tier 1（高优验证组 / High Priority）**：**CHEMBL6001037 (pose 2)**。该分子在 100 ns MD 模拟中展现了口袋内的柔性适应与稳定结合（未解离，RMSD ~0.4–0.65 nm），MM-PBSA 结合自由能极强（ΔG ≈ -37 kcal/mol），且具备优秀的肠道吸收（Caco-2）。建议作为**优先推进湿实验**的候选分子。
    **CHEMBL6001037 (pose 2)** is ranked as Tier 1 due to its stable binding (RMSD ~0.4–0.65 nm without dissociation), strong MM-PBSA binding free energy (ΔG ≈ -37 kcal/mol), and excellent Caco-2 absorption.
*   **Tier 2（骨架参考组 / Scaffold Reference）**：**CHEMBL5571836**。虽然其对接打分全场最高（-11.57 kcal/mol），且 MD 中刚性强、稳定性好，但因其 **QSAR 预测活性低、hERG 心脏毒性极高（0.90）以及 CYP3A4 代谢过快（0.88）**，不建议作为直接成药分子。建议保留其核心骨架作为高亲和力药效团参考。
    **CHEMBL5571836** is classified as Tier 2 (scaffold reference) due to low QSAR prediction, high hERG risk (0.90), and fast metabolism, despite its top docking score and rigid MD stability.
*   **Tier 3（淘汰/避坑组 / Deprioritized）**：结合力较强但 ADMET 极差（如半衰期过短、代谢过快）的分子，建议直接剔除。

### 2. 分子优化方案 / Molecular Optimization Strategy
针对排名前列的化合物普遍存在的高 hERG 心脏毒性及代谢问题，结合构效关系（SAR）提出以下改造建议：
*   **降低 hERG 毒性（多参数优化 MPO）**：hERG 毒性通常与分子的高脂溶性（LogP）及碱性胺基团相关。建议在保留核心骨架的基础上，采用**生物电子等排（Bioisosterism）**替换，引入亲水性基团（如吗啉、哌嗪）以降低 pKa 和脂溶性。
*   **改善代谢稳定性**：针对半衰期较短的分子，可通过封闭易代谢位点（如用氟原子取代易氧化的甲基或芳环），提升体内代谢稳定性。
    To mitigate hERG cardiotoxicity and poor metabolism, we suggest introducing hydrophilic groups (e.g., morpholine) via bioisosterism and blocking metabolic hotspots (e.g., fluorination).

## 🔭 后续优化方向 / Future Directions
1.  **更精确的自由能计算 (FEP/TI)**：当前采用 MM-PBSA 进行快速热力学评估。对于后续确定的同系物先导化合物系列，将引入 **OpenFE** 或 **BindFlow** 等开源工具进行 FEP/TI 计算，以获得更精确的 ΔΔG。
2.  **生成式分子优化**：以本项目筛选出的高潜力分子为起点，利用 **REINVENT 4** 等强化学习生成模型，结合本项目构建的 **QSAR/ADMET 模型作为多参数奖励函数**，进行骨架装饰与 Linker 设计，实现“生成-评估-优化”的闭环。
3.  **蛋白结构精修**：在具备商业许可的前提下，可引入 **Rosetta** 对目标蛋白进行全原子结构精修与突变体建模，进一步提升初始结构质量。

## 🚀 快速运行指南 / Quick Start
1.  安装依赖环境 / Install dependencies: `pip install -r requirements.txt`
2.  运行数据清洗与配体准备 / Prepare ligands: `python 01_data_prep.py`
3.  运行高通量分子对接 / Run docking: `python 02_batch_docking.py`
4.  训练 QSAR 模型 / Train QSAR model: `python 03_qsar_model.py`
5.  分析 MD 轨迹并出图 / Analyze MD: `python 04_md_analysis.py`
6.  ADMET 预测与 MM-PBSA 说明 / ADMET & MM-PBSA: `python 05_admet_mmpbsa.py`

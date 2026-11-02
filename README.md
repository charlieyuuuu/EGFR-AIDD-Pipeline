# EGFR 靶点新型抑制剂的 AI 辅助虚拟筛选与分子动力学验证
# AI-Assisted Virtual Screening & MD Validation Pipeline for EGFR Inhibitors

## 📌 项目概述 / Project Overview
本项目展示了一个完整的 AIDD/CADD（人工智能与计算机辅助药物设计）流水线，用于发现靶向 EGFR 的新型小分子抑制剂。工作流涵盖了**基于结构的虚拟筛选**、**机器学习（QSAR）交叉验证**以及 **100 纳秒（ns）全原子分子动力学（MD）模拟**。核心目的是剔除传统分子对接的假阳性结果，并评估候选复合物的动态结合稳定性。
This repository presents an end-to-end AIDD/CADD pipeline for discovering novel EGFR inhibitors. The workflow integrates **structure-based virtual screening**, **machine learning (QSAR) cross-validation**, and **100 ns all-atom molecular dynamics (MD) simulations**. The primary goal is to eliminate false positives from traditional docking and evaluate the dynamic binding stability.

## 🛠️ 技术栈 / Tech Stack
*   **机器学习与化学信息学 (AIDD) / AI & Cheminformatics:** Scikit-learn (随机森林/Random Forest), RDKit, Pandas (Morgan 分子指纹/Fingerprints)
*   **结构生物学与虚拟筛选 (CADD) / Structure-Based Design:** AutoDock Vina (v1.2.5), Meeko, OpenBabel
*   **分子动力学模拟 (MD) / Molecular Dynamics:** GROMACS 2025 (GPU 加速/accelerated), acpype (GAFF2 力场/force field)
*   **自动化与可视化 / Automation & Visualization:** Linux Shell, Python `subprocess`, Matplotlib, PyMOL

## 📂 研发流水线概述 / Pipeline Overview
1.  **数据准备 / Data Preparation**：从 ChEMBL 数据库中提取 10,734 个具有 EGFR 活性的化合物。利用 RDKit 进行 Lipinski 类药性五规则过滤，并生成 3D 构象及 PDBQT 格式文件。
    Filtered 10,734 EGFR-active compounds from the ChEMBL database using Lipinski's Rule of Five via RDKit. Generated 3D conformers and PDBQT formats.
2.  **高通量虚拟筛选 (CADD) / High-Throughput Virtual Screening**：在 Linux 服务器上编写 Python 自动化脚本，调用 AutoDock Vina 对约 7,500 个化合物完成批量化对接，并根据结合能初筛出 Top 50 候选分子。
    Automated the docking of ~7,500 compounds against the EGFR kinase domain using Python `subprocess` scripts. Selected the Top 50 candidates based on binding affinity.
3.  **AI / QSAR 交叉验证 (AIDD) / AI Cross-Validation**：提取 2048 位 ECFP4 Morgan 分子指纹，训练随机森林回归模型（测试集 R² = 0.682）预测 pIC50。利用该机器学习模型对 Top 50 对接结果进行活性交叉验证，识别假阳性。
    Extracted 2048-bit ECFP4 Morgan fingerprints and trained a Random Forest Regressor (Test R² = 0.682) to predict bioactivity (pIC50). This ML model was used to cross-validate the Top 50 docking hits.
4.  **MD 体系搭建与模拟 (CADD) / MD Setup & Simulation**：针对 6 个代表性候选分子，利用 GROMACS 开展 100 ns 全原子动力学模拟。编写自动化脚本完成拓扑合并（处理 GAFF2 力场参数），并使用 GPU 大幅加速计算。
    Conducted 100 ns MD simulations on 6 representative complexes using GROMACS. Developed automated scripts for topology merging and utilized GPU acceleration.
5.  **轨迹分析 / Trajectory Analysis**：利用 `gmx trjconv`、`gmx rms` 和 `gmx mindist` 计算配体 RMSD、蛋白骨架 RMSD 以及蛋白-配体最小距离，评估复合物的动态稳定性。
    Evaluated complex stability by calculating ligand RMSD, protein backbone RMSD, and minimum distances using GROMACS analysis tools.

## 🔬 核心科学发现 / Key Findings
*静态对接和单一机器学习预测均存在局限性。通过结合 AI 模型与基于物理原理的 MD 模拟，我们揭示了经典的“假阳性”案例：*
*Static docking and pure ML predictions have limitations. By combining AI and physics-based methods, we revealed a critical "false positive" case:*
*   Vina 打分最高（-11.57 kcal/mol）的 **CHEMBL6001037** 主构象（pose 1），在 MD 模拟进行到 3 ns 时即从结合口袋中**解离**。
    The top-ranked docking pose of **CHEMBL6001037** (Vina score: -11.57 kcal/mol) **dissociated** from the binding pocket within 3 ns during the MD simulation.
*   然而，通过多构象筛选出的次要构象（pose 2）在 100 ns 的轨迹中保持了**极高的稳定性**（配体 RMSD 稳定在 0.4 nm 左右）。
    However, a secondary docking pose (pose 2) remained **highly stable** throughout the 100 ns trajectory (ligand RMSD ~ 0.4 nm).
*   *结论：这证明了在早期药物发现中，静态打分绝不能替代动态验证。将 QSAR 过滤与 MD 模拟深度结合，是现代 AIDD 工作流中确保苗头化合物（Hit）可靠性的关键。*
    *Conclusion: This highlights the necessity of dynamic validation over static docking scores. The integration of QSAR filtering and MD simulations is essential in modern AIDD workflows.*

## 📊 结果展示 / Results
![Ligand RMSD](results/Fig1_Ligand_RMSD.png)
*图 1：6 个候选复合物在 100 ns MD 模拟过程中的配体 RMSD 变化曲线。*
*Figure 1: Ligand RMSD over 100 ns MD simulations for 6 candidate complexes.*

## 🚀 快速运行指南 / Quick Start
1.  安装依赖环境 / Install dependencies: `pip install -r requirements.txt`
2.  运行数据清洗与配体准备 / Prepare ligands: `python 01_data_prep.py`
3.  运行高通量分子对接 / Run docking: `python 02_batch_docking.py`
4.  训练 QSAR 模型 / Train QSAR model: `python 03_qsar_model.py`
5.  分析 MD 轨迹并出图 / Analyze MD: `python 04_md_analysis.py`

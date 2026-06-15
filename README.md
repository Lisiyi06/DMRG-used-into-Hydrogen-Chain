# Simulating Electronic Structure of Hydrogen Chains Using Density Matrix Renormalization Group

**课题名称**  
Simulating Ground States of Small Quantum Many-Body Systems Using Classical Numerical Methods  
—— *Simulating Electronic Structure of Hydrogen Chains Using Density Matrix Renormalization Group*

---

## 1. 研究思路

量子多体系统的基态计算是凝聚态物理和量子化学的核心问题。一维强关联系统因量子涨落显著，往往没有解析解，必须借助数值方法。在众多经典算法中，密度矩阵重整化群（DMRG）因其对一维及准一维系统基态描述的高精度，具有显著优势。

本研究从以下三个层面展开：

1. **从零实现 DMRG**：用 Python 实现一个简约但功能完整的 DMRG 程序，深入理解其变分扫荡机制；
2. **横向对比**：在同一物理系统（一维海森堡链）上，对比 DMRG 与 TEBD / iTEBD 的精度、效率与收敛行为；
3. **应用迁移**：将 DMRG 应用于量子化学，计算 H₂ 及线性氢链 Hₙ 的基态能量与键级，从物理层面解释共价键起源，并探讨氢链是否存在稳定的金属相。

---

## 2. 研究内容与技术路线

### 2.1 算法原理与 Python 实现（手写 DMRG）

- 哈密顿量以**矩阵乘积算符（MPO）**形式自动生成，支持海森堡模型、Hubbard 模型及量子化学哈密顿量；
- 波函数表示为**矩阵乘积态（MPS）**，实现左右正则化（混合正则形式）；
- 构造局部二体有效哈密顿量，并利用 **Lanczos 对角化**求解基态；
- 通过**奇异值分解（SVD）**进行截断，并自适应控制键维度。

> **注**：部分高性能 DMRG 包的底层使用 Julia 实现，但为了清晰演示手动复现环节，本项目统一使用 Python。尽管计算速度有一定劣势，但更便于理解算法细节。

### 2.2 方法对比：DMRG vs. TEBD / iTEBD

**测试系统**：一维反铁磁海森堡自旋链

**对比维度**：
- **基态能量精度**：以精确对角化（ED）为基准，比较 DMRG（扫荡收敛）与虚时间演化（iTEBD）的最终能量误差；
- **收敛速度**：记录 DMRG 扫荡次数与 TEBD / iTEBD 虚时间步数的计算时间；
- **纠缠熵增长处理能力**：对于实时间演化（TEBD）模拟量子 quench 动力学，DMRG 无法直接处理。基态问题 DMRG 占优，动力学问题 TEBD 占优。

### 2.3 应用实例：氢分子与氢链的基态能量与结构计算

#### 2.3.1 H₂ 分子：基准测试与基组效应

- 使用 STO‑3G 最小基组，DMRG 计算 H₂ 基态能量为 **-1.137283834489 Hartree**，与 PySCF FCI 精确值完全一致（误差 < 1e-15），验证了 DMRG 实现的正确性。
- 换用 cc‑pVTZ 基组（含极化函数）后，能量降至 **-1.156572416466 Hartree**，更接近实验值/完全基组极限（约 -1.1745 Hartree），说明剩余误差主要来自基组截断而非 DMRG 算法。

#### 2.3.2 等间距一维氢链（Hₙ，n = 2,4,6,8,10,12,14）

- 对每个 n，在 0.6–1.6 Å 范围内用三分法扫描最优键长，并计算相应的基态能量。
- **热力学稳定性**：所有氢链的最优能量均低于 n 个独立氢原子的总能量，表明等间距氢链在热力学上是亚稳的。
- **键长趋势**：最优键长随原子数增加而增大（从 H₂ 的约 0.74 Å 增至 H₁₄ 的约 1.2 Å），说明链增长会削弱化学键。
- **异常变化**：n=2 → 4 时能量和键长变化尤为剧烈，提示等间距假设可能不是真实基态结构。

#### 2.3.3 H₁₀ 的 Peierls 二聚化：一维氢链的“解体”

- **物理动机**：根据 Peierls 定理，一维半满金属在零温下必然发生二聚化（晶格形变打开能隙，降低电子能量），等间距结构只是亚稳态。
- **计算方案**：允许 H₁₀ 中短键（Rₛ）和长键（Rₗ）交替排列，进行二维网格扫描（Rₛ: 0.6–1.2 Å，Rₗ: 0.8–1.4 Å），寻找能量最低的 (Rₛ, Rₗ) 组合。
- **最优二聚化结构**：
  - 短键 Rₛ* ≈ **0.75 Å**，长键 Rₗ* ≈ **1.50 Å**
  - 二聚化参数 δ = (Rₗ − Rₛ)/R̄ ≈ **0.67**
- **能量对比**：
  - H₁₀ 二聚化最优能量：**-5.6879675956 Hartree**
  - 5 个独立 H₂ 分子（R = 0.74 Å）总能量：**-5.6864191724 Hartree**
  - 差值仅 **−0.00155 Hartree（约 −0.042 eV）**，相对误差 **0.027%**
- **结论**：最优二聚化 H₁₀ 链在能量上几乎等同于 5 个 H₂ 分子，表明长氢链在零温下会自发“解体”为独立的氢分子单元。一维氢链的金属相不可能稳定存在。

---

## 3. 项目结构（拟）

```
DMRG-used-into-Hydrogen-Chain/
│
├── project_full_v1.ipynb # 所有代码集合（主程序）
│
├── environment_full.yml # 完整 Conda 环境配置
├── environment_history.yml # 环境历史记录（备用）
│
├── final_paper # 结题论文
├── midterm_report # 中期报告
├── oral_presentation # 答辩PPT
│
├── STEP1-Junjie Xu/ # 阶段1：DMRG 原理实现
├── STEP2-Jingyi Wang/ # 阶段2：DMRG vs TEBD 对比
├── STEP3-Jiangteng Wang & Siyi Li/ # 阶段3：DMRG计算氢链应用
│
├── dmrg.png # DMRG 示意图
├── dmrg_svd.png # SVD 截断示意图
│
└── README.md # 本文件
```

---

## 4. 运行

```bash
# 创建完整环境（推荐）
conda env create -f environment_full.yml

# 或使用历史版本环境
conda env create -f environment_history.yml
```

---

## 5. 依赖环境

- Python 3.8+
- NumPy / SciPy
- Matplotlib
- PySCF（分子积分）
- block2（DMRG 后端）

---

欢迎联系作者
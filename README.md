# Wood Dielectric Inversion

木材纹理方向对微波介电特性影响的反演研究。基于 Keysight VNA 测得的 S21 透射数据，通过 Fabry-Perot 传输模型和多种全局优化算法，提取木材在 18–44 GHz 频段的复介电常数。

## 目录结构

```
D:\woods\
├── run_wood_inversion.m     # 主反演脚本 (含全部算法实现)
├── load_measurement.m       # Keysight VNA CSV 数据读取器
├── plot_wood_comparison.m   # 对比图生成 (独立运行)
├── n_scan_diagnosis.m       # 频率指数 n 扫描诊断工具 (独立运行)
├── lzmwoods\                # 原始测量数据
│   ├── air1.csv / air2.csv / air3.csv   # 空气透射参考 (用于归一化)
│   ├── 0_垂直.csv / 0_平行.csv          # 大木方 (截面大, 纹理方向对比)
│   ├── 1_垂直.csv / 1_平行.csv          # 小木块 (~120mm 厚)
│   └── results_Wood\                    # 反演结果输出目录
│       ├── Data_*_Proposed.mat          # DE-CMA-ES 混合算法结果
│       ├── Data_*_DE.mat / Data_*_CMA.mat ...
│       └── *.png                        # 诊断图和对比图
└── .gitignore
```

## 数据格式

CSV 文件来自 Keysight VNA，`load_measurement.m` 自动解析以下列:

| 列 | 说明 |
|----|------|
| Freq_Hz_ / Freq(Hz) / Frequency | 频率 (Hz) |
| S21(DB) / S21_DB | S21 幅度 (dB) |
| S21(DEG) / S21_DEG / Phase | S21 相位 (deg) |

所有 CSV 必须有 `BEGIN` 行或直接以表头行开始。空气参考文件命名需匹配 `air*.csv` 模式，材料文件为其余 `.csv`。

## 物理模型

### 损耗正切模型 (木材专用)

$$ε^*(f) = ε_r · (1 - j · \tanδ(f))$$

$$\tanδ(f) = \tanδ_{ref} · (f / f_{ref})^n$$

其中 $f_{ref} = 30$ GHz。优化参数为 **4 维**:

| 参数 | 含义 | 典型范围 |
|------|------|----------|
| ε_r | 相对介电常数实部 | 1.5 – 8.0 |
| log10(tanδ_ref) | 参考频率处损耗正切对数 | -3.0 – -0.1 |
| n | 损耗频率色散指数 | -1.5 – 1.5 |
| d | 样品厚度 (m) | 0.10 – 0.50 |

### Fabry-Perot 传输公式

$$T(f) = E · \frac{1 - R^2}{1 - R^2 · E^2}$$

$$E = e^{-γ · d}, \quad γ = j\frac{ω}{c_0}\sqrt{ε^*}, \quad R = \frac{\sqrt{ε^*} - 1}{\sqrt{ε^*} + 1}$$

### 与旧石材模型的区别

旧模型 `eps* = ε_r - j·σ·f^dex/(ε₀·ω)` 适用于导电/半导体材料。木材是介质材料，用损耗正切模型更合适——损耗来自分子弛豫而非电导率。

## 优化算法 (5 种)

| 算法 | 函数名 | 特点 |
|------|--------|------|
| **DE-CMA-ES (Proposed)** | `Hybrid_DE_CMA_ES_Paper_Solver` | 混合差分进化 + 协方差矩阵自适应，每60代交换精英解 |
| **DE** | `Standard_DE_Solver` | 标准差分进化 |
| **CMA-ES** | `Standard_CMA_ES_Solver` | 协方差矩阵自适应进化策略 |
| **GA** | `Standard_GA_Solver` | 遗传算法 (SBX交叉) |
| **PSO** | `Standard_PSO_Solver` | 粒子群优化 |

所有算法共享 LHS-OBL (拉丁超立方采样 + 对立学习) 初始种群策略。

## 代价函数

$$J = 0.5 · mean(L_{mag}) + 0.25 · mean(L_{phase}) + 0.3 · d_{err}^2 + 0.05 · d_{norm}^2$$

- **Huber 损失**: 小残差用 MSE (L2)，大残差用 MAE (L1)，对异常频点鲁棒
- **群延迟约束**: 通过 IFFT 首达峰估计光学路径 `τ_group = n'·d/c₀`，打破 d 与 tanδ 的简并
- **厚度正则化**: 将厚度约束在用户给定范围内

## 预处理流程

```
原始 S21 ──→ 插值到空气参考频率轴 ──→ S21_mat / S21_air (归一化)
                                    │
                                    ├── IFFT 群延迟估计 (仅用于相位校正, 不做门控!)
                                    │
                                    ├── 截取 |S21| > -80 dB 频段
                                    │
                                    ├── Savitzky-Golay 平滑 (2阶, 窗宽 ~N/8)
                                    │
                                    └── 降采样到 ~1000 频点 ──→ 优化器
```

**重要**: 木材是高损耗材料 (S21 < -60 dB)，时域门控会破坏信号。此脚本不做门控，仅用群延迟校正 2π 相位模糊。

## 使用方法

### 1. 全自动批量反演

```matlab
% 在 run_wood_inversion.m 中设置 WORK_MODE = 1
% 确认数据路径:
%   dataDir = 'D:\woods\lzmwoods';
% 确认厚度边界 (大木方 150-500mm, 小木块 100-140mm):
%   lb_default / ub_default, lb_small / ub_small

run_wood_inversion
```

对每个样品依次运行 5 种算法，结果保存为 `results_Wood/Data_<样品>_<算法>.mat`。

### 2. 生成对比图

```matlab
% 前提: 已运行反演, results_Wood/ 下有 .mat 结果文件
plot_wood_comparison
```

输出 4 张图:
- `Wood_Fitting_Comparison.png` — 实测 vs 拟合 (2×2)
- `Wood_Texture_Comparison.png` — 垂直/平行纹理叠加对比
- `Wood_Convergence.png` — 5 算法收敛曲线
- `Wood_Parameter_BarChart.png` — 介电参数条形图

### 3. n 指数扫描诊断

```matlab
n_scan_diagnosis
```

固定 n 值 (-1.5 到 0.5)，对 ε_r、tanδ、d 做 fmincon 优化，绘制 cost vs n 曲线。用于检查损耗模型的频率依赖性是否合理。

### 4. 仅出汇总表 (不重新跑反演)

```matlab
% 设置 WORK_MODE = 2, 然后:
run_wood_inversion
```

## 依赖

- MATLAB R2018b 或更高版本
- 需以下工具箱:
  - Statistics and Machine Learning Toolbox (`lhsdesign`)
  - Signal Processing Toolbox (`sgolayfilt`)
  - Optimization Toolbox (`fmincon`, 仅 `n_scan_diagnosis.m` 需要)

## 结果解读

每个 `.mat` 文件包含:

| 字段 | 说明 |
|------|------|
| `fGHz` | 频率轴 (GHz) |
| `ym_raw` | 实测 |S21| 幅度 |
| `fit_mag` | 模型拟合 |T(f)| |
| `loss_curve` | 优化过程 cost 下降曲线 |
| `params` | [ε_r, log10(tanδ_ref), n, d] |
| `rmse` | 最终 RMSE (dB) |

**已知局限**: 损耗正切幂律模型当 n ≈ -1 时，衰减系数 α 近似与频率无关，导致 |T(f)| 在频段内近似平坦。实测数据通常有 10+ dB 的频变，当前模型无法完全复现。进一步改进可考虑 **Debye 弛豫模型** (增加弛豫频率参数) 来捕捉木材水分含量的色散特征。

# 木材分层非均匀介质的电磁波传播：完整数学推导

## 1. 几何构型与各向异性本构关系

### 1.1 坐标定义

将木材样品置于自由空间中，定义坐标系：

- **z 轴**：电磁波传播方向（正入射，垂直于木块表面）
- **x-y 平面**：木块截面
- **L 轴**：木材纵向（纤维方向）
- **R 轴**：径向（年轮生长方向）
- **T 轴**：弦向

对于径切面入射的样品（z ∥ R）：

- x ∥ L（纤维方向，面内）
- y ∥ T（弦向，面内）
- z ∥ R（传播方向）

### 1.2 正交各向异性本构关系

频域中（采用工程时谐约定 $e^{j\omega t}$），本构关系为：

$$\mathbf{D}(\mathbf{r}, \omega) = \varepsilon_0 \, \bar{\bar{\varepsilon}}_r(\mathbf{r}, \omega) \cdot \mathbf{E}(\mathbf{r}, \omega) \tag{1.1}$$

$$\bar{\bar{\varepsilon}}_r = \begin{bmatrix} \varepsilon_L^* & 0 & 0 \\ 0 & \varepsilon_R^* & 0 \\ 0 & 0 & \varepsilon_T^* \end{bmatrix} \tag{1.2}$$

其中各主轴分量为复数：

$$\varepsilon_i^* = \varepsilon_i' - j\varepsilon_i'' \quad (i = L, R, T) \tag{1.3}$$

- 实部 $\varepsilon_i'$：表征材料的储能能力
- 虚部 $\varepsilon_i''$：表征介质损耗（电磁波能量转化为热能）

木材为非磁性介质，磁本构关系简化为 $\mathbf{B} = \mu_0 \mathbf{H}$。

### 1.3 准周期分层结构

沿 z 轴（传播方向），木材由早材（Earlywood, E）和晚材（Latewood, L）交替层叠构成：

$$\varepsilon_r^*(z) = \begin{cases} \varepsilon_E^* & n = 1, 3, 5, \ldots \quad \text{（早材层）} \\ \varepsilon_L^* & n = 2, 4, 6, \ldots \quad \text{（晚材层）} \end{cases} \tag{1.4}$$

设总层数为 $N$（取偶数，保证完整的 E-L 层对），则层对数 $P = N/2$。第 $n$ 层厚度为 $d_n$，总厚度约束：

$$\sum_{n=1}^{N} d_n = \frac{N}{2}(d_E + d_L) = D_{\text{total}} \tag{1.5}$$

---

## 2. Maxwell 方程组与波动方程

### 2.1 频域 Maxwell 旋度方程

时谐场 $e^{j\omega t}$ 下，无源区域的 Maxwell 旋度方程组为：

$$\nabla \times \mathbf{E} = -j\omega\mu_0 \mathbf{H} \tag{2.1}$$

$$\nabla \times \mathbf{H} = j\omega\varepsilon_0 \bar{\bar{\varepsilon}}_r \cdot \mathbf{E} \tag{2.2}$$

### 2.2 正入射条件下的退化

正入射时 $\mathbf{k} = k_z \hat{\mathbf{z}}$，平面波解的形式为：

$$\mathbf{E}(\mathbf{r}) = \mathbf{E}_0(z) \, e^{-j k_z z}$$

由 $\partial/\partial x = \partial/\partial y = 0$，从 (2.1) 和 (2.2) 可证 $E_z = 0$、$H_z = 0$。电场和磁场均位于 x-y 平面内（相对于 z 方向的 TEM 波）。

### 2.3 两种特征偏振模式

正入射线偏振波的电场方向决定激发哪个介电主轴分量。

**模式 1：平行纹理（$\parallel$）**——电场沿 L 轴（纤维方向，取为 x 轴）

$$\mathbf{E} = E_x(z) \hat{\mathbf{x}}, \quad \mathbf{H} = H_y(z) \hat{\mathbf{y}} \tag{2.3}$$

由 (2.2) 可知 $H_y$ 仅通过 $\varepsilon_L^*$ 与 $E_x$ 耦合。本构退化为标量：

$$\varepsilon_{\parallel}^* = \varepsilon_L^* \tag{2.4}$$

**模式 2：垂直纹理（$\perp$）**——电场沿 T 轴（垂直于纤维，取为 y 轴）

$$\mathbf{E} = E_y(z) \hat{\mathbf{y}}, \quad \mathbf{H} = -H_x(z) \hat{\mathbf{x}} \tag{2.5}$$

本构退化为标量：

$$\varepsilon_{\perp}^* = \varepsilon_T^* \tag{2.6}$$

### 2.4 关键结论

正入射条件下，正交各向异性张量 $\bar{\bar{\varepsilon}}_r$ 对每种偏振模式退化为一个标量 $\varepsilon^*(z)$。两种模式**相互独立**，可**分别测量、分别反演**。以下统一处理标量问题。

---

## 3. 单层均匀介质中的波传播

### 3.1 标量 Helmholtz 方程

以模式 1（$\parallel$）为例。从 (2.1) 和 (2.2) 消去 $H_y$：

$$\frac{\partial E_x}{\partial z} = -j\omega\mu_0 H_y \tag{3.1a}$$

$$\frac{\partial H_y}{\partial z} = -j\omega\varepsilon_0 \varepsilon_r^* E_x \tag{3.1b}$$

交叉求导得 Helmholtz 方程：

$$\frac{\partial^2 E_x}{\partial z^2} = -j\omega\mu_0 \frac{\partial H_y}{\partial z} = -j\omega\mu_0 \cdot (-j\omega\varepsilon_0\varepsilon_r^*) E_x$$

$$\boxed{\frac{d^2 E_x(z)}{dz^2} + k_0^2 \varepsilon_r^*(z) E_x(z) = 0} \tag{3.2}$$

其中 $k_0 = \dfrac{\omega}{c_0} = \dfrac{2\pi f}{c_0}$ 为自由空间波数，$c_0 = 1/\sqrt{\mu_0\varepsilon_0} \approx 2.998 \times 10^8$ m/s。

### 3.2 通解：前向波与后向波

在第 $n$ 层内，$\varepsilon_r^* = \varepsilon_n^*$ 为常数。(3.2) 是常系数二阶微分方程，其通解为：

$$E_x(z) = A_n e^{-j k_0 n_n^* (z - z_{n-1})} + B_n e^{+j k_0 n_n^* (z - z_{n-1})} \tag{3.3}$$

其中 $n_n^* = \sqrt{\varepsilon_n^*}$ 为第 $n$ 层的**复折射率**，$z_{n-1}$ 为第 $n$ 层前表面坐标。

- **第一项** $A_n e^{-j k_0 n_n^* (z - z_{n-1})}$：沿 $+z$ 方向传播的前向波（入射/透射方向）
- **第二项** $B_n e^{+j k_0 n_n^* (z - z_{n-1})}$：沿 $-z$ 方向传播的后向波（反射波）

复折射率的虚部引起**指数衰减**：

$$n_n^* = n_n' - j n_n'' = \sqrt{\varepsilon_n' - j\varepsilon_n''} \tag{3.4}$$

$$e^{-j k_0 (n_n' - j n_n'') d_n} = e^{-j k_0 n_n' d_n} \cdot \underbrace{e^{-k_0 n_n'' d_n}}_{\text{衰减因子 } < 1} \tag{3.5}$$

### 3.3 磁场表达式

由 (3.1a)：

$$H_y(z) = -\frac{1}{j\omega\mu_0} \frac{\partial E_x}{\partial z}$$

代入 (3.3)：

$$H_y(z) = \frac{n_n^*}{Z_0} \left[ A_n e^{-j k_0 n_n^* (z - z_{n-1})} - B_n e^{+j k_0 n_n^* (z - z_{n-1})} \right] \tag{3.6}$$

其中 $Z_0 = \sqrt{\mu_0/\varepsilon_0} \approx 377\,\Omega$ 为自由空间波阻抗。

定义第 $n$ 层的**介质波阻抗**：

$$\boxed{Z_n = \frac{Z_0}{n_n^*} = \frac{Z_0}{\sqrt{\varepsilon_n^*}}} \tag{3.7}$$

则 (3.6) 重写为：

$$H_y(z) = \frac{1}{Z_n} \left[ A_n e^{-j\phi_n(z)} - B_n e^{+j\phi_n(z)} \right] \tag{3.8}$$

其中 $\phi_n(z) = k_0 n_n^* (z - z_{n-1})$ 为局部相位。

---

## 4. 界面边界条件与单层特征矩阵

### 4.1 界面连续性条件

在任意无源界面 $z = z_n$ 处，Maxwell 方程的积分形式要求：

$$\boxed{E_x(z_n^-) = E_x(z_n^+)} \quad \text{（切向电场连续）} \tag{4.1a}$$

$$\boxed{H_y(z_n^-) = H_y(z_n^+) } \quad \text{（切向磁场连续，无表面电流）} \tag{4.1b}$$

### 4.2 单层特征矩阵的构造

考虑第 $n$ 层，厚度 $d_n$，前表面 $z = z_{n-1}$，后表面 $z = z_n = z_{n-1} + d_n$。

在 $z = z_{n-1}^+$（刚进入第 $n$ 层）：

$$E_x(z_{n-1}^+) = A_n + B_n \tag{4.2}$$

$$H_y(z_{n-1}^+) = \frac{1}{Z_n}(A_n - B_n) \tag{4.3}$$

在 $z = z_n^-$（即将离开第 $n$ 层）：

$$E_x(z_n^-) = A_n e^{-j\delta_n} + B_n e^{+j\delta_n} \tag{4.4}$$

$$H_y(z_n^-) = \frac{1}{Z_n}\left(A_n e^{-j\delta_n} - B_n e^{+j\delta_n}\right) \tag{4.5}$$

其中定义了**复电长度**（复相位延迟）：

$$\boxed{\delta_n = k_0 n_n^* d_n = \frac{2\pi f}{c_0} \sqrt{\varepsilon_n^*} \, d_n} \tag{4.6}$$

由 (4.2)-(4.3) 解出 $A_n, B_n$：

$$A_n = \frac{1}{2}\left[E_x(z_{n-1}^+) + Z_n H_y(z_{n-1}^+)\right] \tag{4.7}$$

$$B_n = \frac{1}{2}\left[E_x(z_{n-1}^+) - Z_n H_y(z_{n-1}^+)\right] \tag{4.8}$$

代入 (4.4)-(4.5)，利用 Euler 公式 $\cos\delta_n = (e^{j\delta_n} + e^{-j\delta_n})/2$ 和 $\sin\delta_n = (e^{j\delta_n} - e^{-j\delta_n})/(2j)$ 整理：

$$E_x(z_{n-1}^+) = E_x(z_n^-) \cos\delta_n + j Z_n H_y(z_n^-) \sin\delta_n \tag{4.9}$$

$$H_y(z_{n-1}^+) = \frac{j}{Z_n} E_x(z_n^-) \sin\delta_n + H_y(z_n^-) \cos\delta_n \tag{4.10}$$

写成矩阵形式，定义第 $n$ 层的**特征矩阵**（从后表面映射到前表面）：

$$\boxed{\begin{bmatrix} E_x(z_{n-1}^+) \\ H_y(z_{n-1}^+) \end{bmatrix} = \mathbf{M}_n \begin{bmatrix} E_x(z_n^-) \\ H_y(z_n^-) \end{bmatrix}} \tag{4.11}$$

$$\boxed{\mathbf{M}_n = \begin{bmatrix} \cos\delta_n & j Z_n \sin\delta_n \\ \dfrac{j}{Z_n} \sin\delta_n & \cos\delta_n \end{bmatrix}} \tag{4.12}$$

### 4.3 行列式恒为 1

$$\det(\mathbf{M}_n) = \cos^2\delta_n - j^2 \sin^2\delta_n = \cos^2\delta_n + \sin^2\delta_n = 1 \tag{4.13}$$

每层特征矩阵的行列式恒为 1，与频率、介电常数、厚度均无关。这一性质保证了后续 Chebyshev 加速的数学合法性，也保证了数值稳定性。

---

## 5. 多层介质的总转移矩阵

### 5.1 级联乘积

层间界面处电磁场连续，各层矩阵按空间顺序直接连乘。输入端（$z = 0$，第 1 层前表面）到输出端（$z = D_{\text{total}}$，第 $N$ 层后表面）：

$$\begin{bmatrix} E_x(0) \\ H_y(0) \end{bmatrix} = \mathbf{M}_1 \mathbf{M}_2 \cdots \mathbf{M}_N \begin{bmatrix} E_x(D_{\text{total}}) \\ H_y(D_{\text{total}}) \end{bmatrix} \tag{5.1}$$

定义**总转移矩阵**：

$$\boxed{\mathbf{M}_{\text{total}} = \prod_{n=1}^{N} \mathbf{M}_n = \begin{bmatrix} A(f) & B(f) \\ C(f) & D(f) \end{bmatrix}} \tag{5.2}$$

对于早材/晚材交替结构：

$$\mathbf{M}_n = \begin{cases} \mathbf{M}_E = \begin{bmatrix} \cos\delta_E & j Z_E \sin\delta_E \\ \dfrac{j}{Z_E}\sin\delta_E & \cos\delta_E \end{bmatrix} & n \text{ 为奇数} \\[12pt] \mathbf{M}_L = \begin{bmatrix} \cos\delta_L & j Z_L \sin\delta_L \\ \dfrac{j}{Z_L}\sin\delta_L & \cos\delta_L \end{bmatrix} & n \text{ 为偶数} \end{cases} \tag{5.3}$$

其中：

$$\delta_E = k_0 \sqrt{\varepsilon_E^*} \, d_E, \quad Z_E = \frac{Z_0}{\sqrt{\varepsilon_E^*}} \tag{5.4a}$$

$$\delta_L = k_0 \sqrt{\varepsilon_L^*} \, d_L, \quad Z_L = \frac{Z_0}{\sqrt{\varepsilon_L^*}} \tag{5.4b}$$

### 5.2 物理诠释

每一步矩阵乘法 $\mathbf{M}_n$ 完成两项物理操作：
1. **相位累积**：$\cos\delta_n$ 和 $\sin\delta_n$ 记录了波在第 $n$ 层内的传播相位和衰减
2. **阻抗变换**：$Z_n$ 和 $1/Z_n$ 通过正弦项耦合了电场与磁场，表征了界面处的阻抗不连续性导致的局部反射

总矩阵 $\mathbf{M}_{\text{total}}$ 的四个元素 $A, B, C, D$ 完整编码了所有层间多重反射和透射的相干叠加。

### 5.3 Chebyshev 多项式加速

当层数 $N$ 很大时（例如 120 mm 木块，年轮间距 5 mm，$N = 48$），每个频点做 $N$ 次 $2 \times 2$ 复数矩阵乘法在优化迭代中代价高昂。利用 $\det(\mathbf{M}_{\text{pair}}) = 1$ 的性质，可用 **Chebyshev 恒等式**将计算复杂度从 $O(N)$ 降为 $O(1)$。

定义单层对转移矩阵：

$$\mathbf{M}_{\text{pair}} = \mathbf{M}_E \cdot \mathbf{M}_L = \begin{bmatrix} \alpha & \beta \\ \gamma & \delta \end{bmatrix} \tag{5.5}$$

由各层 $\det = 1$ 及行列式乘积性质，有：

$$\det(\mathbf{M}_{\text{pair}}) = \alpha\delta - \beta\gamma = 1 \tag{5.6}$$

定义半迹（矩阵迹的一半）：

$$s = \frac{\alpha + \delta}{2} = \frac{1}{2} \operatorname{tr}(\mathbf{M}_{\text{pair}}) \tag{5.7}$$

**Chebyshev 恒等式**（可由特征值分解直接证明）：对于任意满足 $\det = 1$ 的 $2 \times 2$ 矩阵，其 $P$ 次幂可用第二类 Chebyshev 多项式表达：

$$\boxed{\mathbf{M}_{\text{pair}}^P = U_{P-1}(s) \cdot \mathbf{M}_{\text{pair}} - U_{P-2}(s) \cdot \mathbf{I}} \tag{5.8}$$

其中 $U_n(s)$ 为第二类 Chebyshev 多项式，满足递推关系 $U_{n+1}(s) = 2s U_n(s) - U_{n-1}(s)$，初始条件 $U_{-1}(s) \equiv 0$，$U_0(s) \equiv 1$。

### 5.4 $U_n(s)$ 的计算

根据 $s$ 的取值范围，采用不同的数值计算策略：

**情况 1：传播模式（$|s| \le 1$）**

令 $s = \cos\theta$，即 $\theta = \arccos(s)$。Chebyshev 多项式的三角形式为：

$$U_n(\cos\theta) = \frac{\sin[(n+1)\theta]}{\sin\theta} \tag{5.9}$$

**情况 2：倏逝模式（$|s| > 1$）**

令 $s = \cosh\xi$，即 $\xi = \operatorname{arccosh}(s) = \ln\left(s + \sqrt{s^2 - 1}\right)$。双曲形式为：

$$U_n(\cosh\xi) = \frac{\sinh[(n+1)\xi]}{\sinh\xi} \tag{5.10}$$

**统一实现**：在 MATLAB 中，$\arccos$ 对任意复数参数均有定义，公式 (5.9) 自动涵盖 (5.10)。因此实现时可直接使用：

```matlab
theta = acos(s);                          % 对任意复数 s 均有效
U_Pm1 = sin(P * theta) ./ sin(theta);     % U_{P-1}(s)
U_Pm2 = sin((P-1) * theta) ./ sin(theta); % U_{P-2}(s)
M_total = U_Pm1 * M_pair - U_Pm2 * eye(2); % Chebyshev 恒等式
```

### 5.5 每频点的计算步骤

1. 计算复折射率与电长度：

$$\delta_E = k_0 \sqrt{\varepsilon_E^*} d_E, \quad \delta_L = k_0 \sqrt{\varepsilon_L^*} d_L \tag{5.11}$$

2. 构建早材和晚材的特征矩阵 $\mathbf{M}_E$, $\mathbf{M}_L$
3. 计算层对矩阵 $\mathbf{M}_{\text{pair}} = \mathbf{M}_E \cdot \mathbf{M}_L$
4. 计算半迹 $s = (\alpha + \delta)/2$
5. $\theta = \arccos(s)$
6. $U_{P-1} = \sin(P\theta)/\sin\theta$，$U_{P-2} = \sin((P-1)\theta)/\sin\theta$
7. $\mathbf{M}_{\text{total}} = U_{P-1} \cdot \mathbf{M}_{\text{pair}} - U_{P-2} \cdot \mathbf{I}$

**计算量与层数 $N$ 无关。**

---

## 6. S 参数前向解析表达式

### 6.1 从转移矩阵到 S 参数

木材两侧为无限大自由空间（波阻抗 $Z_0$）。

在 $z = 0^-$（自由空间侧，输入端口）：

$$E_x(0^-) = E_{\text{inc}} + E_{\text{ref}} \tag{6.1}$$

$$H_y(0^-) = \frac{1}{Z_0}(E_{\text{inc}} - E_{\text{ref}}) \tag{6.2}$$

在 $z = D_{\text{total}}^+$（自由空间侧，输出端口）：

$$E_x(D_{\text{total}}^+) = E_{\text{trans}} \tag{6.3}$$

$$H_y(D_{\text{total}}^+) = \frac{1}{Z_0} E_{\text{trans}} \tag{6.4}$$

界面连续性条件：$E_x(0^-) = E_x(0^+)$ 和 $H_y(0^-) = H_y(0^+)$（在 $z = 0$ 处），以及在 $z = D_{\text{total}}$ 处的对应条件。

将 (6.1)-(6.2) 代入总转移矩阵关系 (5.1) 的左侧，将 (6.3)-(6.4) 代入右侧：

$$\begin{bmatrix} E_{\text{inc}} + E_{\text{ref}} \\ \dfrac{1}{Z_0}(E_{\text{inc}} - E_{\text{ref}}) \end{bmatrix} = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \begin{bmatrix} E_{\text{trans}} \\ \dfrac{1}{Z_0}E_{\text{trans}} \end{bmatrix} \tag{6.5}$$

### 6.2 S21 的推导

展开 (6.5) 的第一行：

$$E_{\text{inc}} + E_{\text{ref}} = \left(A + \frac{B}{Z_0}\right) E_{\text{trans}} \tag{6.6}$$

展开第二行：

$$\frac{1}{Z_0}(E_{\text{inc}} - E_{\text{ref}}) = \left(C + \frac{D}{Z_0}\right) E_{\text{trans}} \tag{6.7}$$

将 (6.6) 与 (6.7) 相加：

$$2E_{\text{inc}} = \left(A + \frac{B}{Z_0} + C Z_0 + D\right) E_{\text{trans}} \tag{6.8}$$

由此得透射系数：

$$\boxed{S_{21}(f) = \frac{E_{\text{trans}}}{E_{\text{inc}}} = \frac{2}{A + \dfrac{B}{Z_0} + C Z_0 + D}} \tag{6.9}$$

### 6.3 S11 的推导

将 (6.6) 与 (6.7) 相减：

$$2E_{\text{ref}} = \left(A + \frac{B}{Z_0} - C Z_0 - D\right) E_{\text{trans}} \tag{6.10}$$

两边除以 $2E_{\text{inc}}$ 并代入 (6.9)：

$$\boxed{S_{11}(f) = \frac{E_{\text{ref}}}{E_{\text{inc}}} = \frac{A + \dfrac{B}{Z_0} - C Z_0 - D}{A + \dfrac{B}{Z_0} + C Z_0 + D}} \tag{6.11}$$

其中 $A(f), B(f), C(f), D(f)$ 为总转移矩阵 $\mathbf{M}_{\text{total}}$ 的频率相关元素。

### 6.4 物理诠释

- **$S_{21}$（透射系数）**：包含幅度衰减和相位延迟。对于 120 mm 厚的湿木材，$|S_{21}|$ 在 18-44 GHz 约为 -20 到 -60 dB。
- **$S_{11}$（反射系数）**：表征各层间阻抗失配引起的总反射。对于木材（$\sqrt{\varepsilon^*} \approx 2\text{-}3$），界面反射较弱（$R \approx -0.2$ 到 $-0.3$）。

当前实验仅有 $S_{21}$ 数据。若同时获取 $S_{11}$，两条曲线联合约束可**成倍增加信息量**，显著改善反演的可辨识性。

---

## 7. 特殊情况验证

### 7.1 均匀介质退化极限

若所有层的介电常数相同（$\varepsilon_E^* = \varepsilon_L^* = \varepsilon^*$），则 $\mathbf{M}_E = \mathbf{M}_L = \mathbf{M}$。由 Chebyshev 恒等式，$\mathbf{M}_{\text{pair}} = \mathbf{M}^2$，$s = \cos 2\delta$：

$$\mathbf{M}^{2P} = \begin{bmatrix} \cos(2P\delta) & j Z \sin(2P\delta) \\ \dfrac{j}{Z}\sin(2P\delta) & \cos(2P\delta) \end{bmatrix} \tag{7.1}$$

代入 (6.9)，经过代数化简即得到标准的 Fabry-Perot 透射公式：

$$S_{21} = \frac{(1-R^2) e^{-j k_0 n^* D_{\text{total}}}}{1 - R^2 e^{-j 2 k_0 n^* D_{\text{total}}}}, \quad R = \frac{Z - Z_0}{Z + Z_0} = \frac{1/\sqrt{\varepsilon^*} - 1}{1/\sqrt{\varepsilon^*} + 1} \tag{7.2}$$

这验证了分层 TMM 推导在均匀极限下的一致性——Fabry-Perot 公式是本推导的特例。

### 7.2 高损耗极限

当 $\varepsilon''$ 很大时（湿木材的典型情况），$\delta_n$ 有较大的虚部，$|\sin\delta_n| \gg 1$ 且 $|\cos\delta_n| \gg 1$。Fabry-Perot 谐振条纹被完全压制。$|S_{21}|$ 主要由指数衰减 $e^{-k_0 \bar{n}'' D_{\text{total}}}$ 主导，频谱呈现平滑衰减曲线。

---

## 8. 逆问题：参数反演

### 8.1 前向映射

完整的前向映射将物理参数空间映射到可测量 S 参数空间：

$$\boxed{S_{21}(f) = \mathcal{F}\left(\varepsilon_E'(f), \varepsilon_E''(f), \varepsilon_L'(f), \varepsilon_L''(f), d_E, d_L, N; \; f\right)} \tag{8.1}$$

$$\boxed{S_{11}(f) = \mathcal{G}\left(\varepsilon_E'(f), \varepsilon_E''(f), \varepsilon_L'(f), \varepsilon_L''(f), d_E, d_L, N; \; f\right)} \tag{8.2}$$

### 8.2 Debye 弛豫频率色散模型

木材中的水分以两种状态存在：自由水（细胞腔内）和结合水（细胞壁内）。每种水分贡献一个 Debye 弛豫过程。对于单弛豫近似：

$$\boxed{\varepsilon^*(f) = \varepsilon_\infty + \frac{\Delta\varepsilon}{1 + j f / f_{\text{relax}}}} \tag{8.3}$$

展开为实部和虚部：

$$\varepsilon'(f) = \varepsilon_\infty + \frac{\Delta\varepsilon}{1 + (f/f_{\text{relax}})^2} \tag{8.4a}$$

$$\varepsilon''(f) = \Delta\varepsilon \cdot \frac{f/f_{\text{relax}}}{1 + (f/f_{\text{relax}})^2} \tag{8.4b}$$

参数物理意义：
- $\varepsilon_\infty$：高频极限介电常数（$f \gg f_{\text{relax}}$ 时的干木材基质）
- $\Delta\varepsilon$：弛豫强度（正比于水分含量）
- $f_{\text{relax}}$：弛豫特征频率（自由水 $\sim$18 GHz，结合水 $\sim$1-10 GHz）

早材和晚材具有**相同的弛豫频率**（水分类型相同），但**不同的 $\varepsilon_\infty$ 和 $\Delta\varepsilon$**（水分含量和密度不同）。

> **物理机制说明**：微观上早材与晚材的微观孔隙结构及结合水/自由水比例存在差异。但在 18–40 GHz 高频段，电磁响应的绝对主导极化机制为自由水的强 Debye 弛豫，结合水在该频段的贡献已基本趋于饱和并归入高频背景介电常数 $\varepsilon_\infty$ 中。因此，本模型将早、晚材的主导弛豫频率统一等效为自由水在木质基质中的平均弛豫频率 $f_{\text{relax}}$。这一近似在保证高频水分响应物理机理的同时，使寻优变量成功减少 1 个，显著提升了逆问题反演的收敛速度与解的唯一性。

### 8.3 优化问题

给定 $M$ 个频率点的测量数据 $S_{21}^{\text{meas}}(f_m)$，反演问题表述为：

$$\boxed{\mathbf{X}^* = \arg\min_{\mathbf{X}} \sum_{m=1}^{M} \left| S_{21}^{\text{meas}}(f_m) - \mathcal{F}(f_m; \mathbf{X}) \right|^2} \tag{8.5}$$

对于 $S_{11}$ 和 $S_{21}$ 联合反演：

$$J(\mathbf{X}) = \sum_{m=1}^{M} \left( w_1 \left|S_{21}^{\text{meas}} - S_{21}^{\text{model}}\right|^2 + w_2 \left|S_{11}^{\text{meas}} - S_{11}^{\text{model}}\right|^2 \right) \tag{8.6}$$

### 8.4 参数降维策略

| 参数 | 降维方案 | 依据 |
|------|---------|------|
| $d_E$, $d_L$ | 固定比值（年轮宽度测量） | 可从样品截面图像精确测量 |
| $N$ | $N = \text{round}(D_{\text{total}} / (d_E + d_L))$ | 由总厚度和层对间距完全确定 |
| $D_{\text{total}}$ | 卡尺测量，固定 | 独立于电磁测量的几何约束 |
| $\Delta\varepsilon_E$, $\Delta\varepsilon_L$ | 可固定（文献值）或加入优化 | 由水分含量的独立测量约束 |
| $\varepsilon_{\infty,E}$, $\varepsilon_{\infty,L}$, $f_{\text{relax}}$ | **优化变量（3 个）** | 核心未知量 |

---

## 9. 核心公式汇总

| 物理量 | 表达式 |
|--------|--------|
| Debye 复介电常数 | $\varepsilon^*(f) = \varepsilon_\infty + \dfrac{\Delta\varepsilon}{1 + j f / f_{\text{relax}}}$ |
| 复折射率 | $n^* = \sqrt{\varepsilon^*}$ |
| 电长度 | $\delta_n = \dfrac{2\pi f}{c_0} \sqrt{\varepsilon_n^*} \, d_n$ |
| 介质波阻抗 | $Z_n = \dfrac{Z_0}{\sqrt{\varepsilon_n^*}}$ |
| 单层特征矩阵 | $\mathbf{M}_n = \begin{bmatrix} \cos\delta_n & j Z_n \sin\delta_n \\ \dfrac{j}{Z_n}\sin\delta_n & \cos\delta_n \end{bmatrix}$ |
| 层对矩阵幂（Chebyshev） | $\mathbf{M}_{\text{pair}}^P = U_{P-1}(s) \mathbf{M}_{\text{pair}} - U_{P-2}(s) \mathbf{I}$ |
| 半迹 | $s = \dfrac{\alpha + \delta}{2}$ |
| Chebyshev 多项式 | $U_n(s) = \dfrac{\sin[(n+1)\arccos(s)]}{\sin[\arccos(s)]}$ |
| 透射系数 | $S_{21} = \dfrac{2}{A + B/Z_0 + C Z_0 + D}$ |
| 反射系数 | $S_{11} = \dfrac{A + B/Z_0 - C Z_0 - D}{A + B/Z_0 + C Z_0 + D}$ |

所有系数 $A, B, C, D$ 均为频率 $f$ 的函数，由总转移矩阵 $\mathbf{M}_{\text{total}}$ 的元素给出。

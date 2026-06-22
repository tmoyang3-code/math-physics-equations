# 第3章 波动方程初值问题与行波法

> ⭐ **考试重点**：达朗贝尔公式为必考解答题（12分），齐次化原理和泊松公式为高频考点

---

## 目录

- [3.1 一维初值问题（达朗贝尔公式）](#31-一维初值问题达朗贝尔公式)
- [3.2 三维初值问题（球面平均法、泊松公式）](#32-三维初值问题球面平均法泊松公式)
- [3.3 二维初值问题（降维法）](#33-二维初值问题降维法)
- [3.4 依赖区域、决定区域、影响区域](#34-依赖区域决定区域影响区域)
- [📝 考试真题与解题模板](#-考试真题与解题模板)

---

## 3.1 一维初值问题（达朗贝尔公式）

### 3.1.1 齐次方程（达朗贝尔公式）

考虑无限长弦的自由振动问题：

$$
\begin{cases}
u_{tt} = a^2 u_{xx}, & x \in \mathbb{R}^1, \; t > 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & x \in \mathbb{R}^1
\end{cases} \tag{1}
$$

**行波法求解：**

作特征变换：
$$
\xi = x + at, \quad \eta = x - at
$$

方程化为：
$$
u_{\xi\eta} = 0 \tag{2}
$$

其通解为：
$$
u(x,t) = f(\xi) + g(\eta) = f(x + at) + g(x - at) \tag{3}
$$

其中 $f, g$ 为任意二次连续可微函数，分别表示**右行波**和**左行波**。

---

**利用初值条件确定 $f, g$：**

$$
u(x,0) = \varphi(x) = f(x) + g(x) \tag{4}
$$

$$
u_t(x,0) = \psi(x) = a f'(x) - a g'(x) \tag{5}
$$

对 (5) 式积分得：
$$
f(x) - g(x) = \frac{1}{a} \int_{x_0}^x \psi(s) ds + C \tag{6}
$$

其中 $x_0, C$ 为任意常数。联立 (4) 和 (6) 解得：

$$
f(x) = \frac{1}{2}\varphi(x) + \frac{1}{2a} \int_{x_0}^x \psi(s) ds + \frac{C}{2} \tag{7}
$$

$$
g(x) = \frac{1}{2}\varphi(x) - \frac{1}{2a} \int_{x_0}^x \psi(s) ds - \frac{C}{2} \tag{8}
$$

---

### ✨ 达朗贝尔公式（必背）

$$
u(x,t) = \frac{1}{2}\left[\varphi(x + at) + \varphi(x - at)\right] + \frac{1}{2a} \int_{x - at}^{x + at} \psi(s) ds \tag{9}
$$

> 💡 **物理意义**：初始位移 $\varphi(x)$ 分为两个波，分别以速度 $a$ 向左右传播；初始速度 $\psi(x)$ 的影响表现为积分形式。

---

### 3.1.2 非齐次方程与齐次化原理

考虑非齐次初值问题：

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + f(x,t), & x \in \mathbb{R}^1, \; t > 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & x \in \mathbb{R}^1
\end{cases} \tag{10}
$$

**叠加原理**：令 $u = v + w$，其中：

- $v$ 满足齐次方程非齐次初值：
$$
\begin{cases}
v_{tt} = a^2 v_{xx}, & x \in \mathbb{R}^1, \; t > 0 \\
v(x,0) = \varphi(x), \; v_t(x,0) = \psi(x), & x \in \mathbb{R}^1
\end{cases} \tag{11}
$$

- $w$ 满足非齐次方程齐次初值：
$$
\begin{cases}
w_{tt} = a^2 w_{xx} + f(x,t), & x \in \mathbb{R}^1, \; t > 0 \\
w(x,0) = 0, \; w_t(x,0) = 0, & x \in \mathbb{R}^1
\end{cases} \tag{12}
$$

---

### 齐次化原理（Duhamel原理）

**定理1**：设 $f(x,t) \in C^1$，若 $h = h(x,t;\tau)$ 是问题

$$
\begin{cases}
h_{tt} = a^2 h_{xx}, & x \in \mathbb{R}^1, \; t > \tau \\
h|_{t = \tau} = 0, \; h_t|_{t = \tau} = f(x,\tau), & x \in \mathbb{R}^1
\end{cases} \tag{15}
$$

的解，则积分

$$
w(x,t) = \int_0^t h(x,t;\tau) d\tau \tag{16}
$$

是问题 (12) 的解。

**推导**：作平移变换 $t' = t - \tau$，应用达朗贝尔公式：

$$
h(x,t;\tau) = \frac{1}{2a} \int_{x - a(t - \tau)}^{x + a(t - \tau)} f(\xi,\tau) d\xi
$$

因此：

$$
w(x,t) = \frac{1}{2a} \int_0^t \int_{x - a(t - \tau)}^{x + a(t - \tau)} f(\xi,\tau) d\xi d\tau
$$

---

### 非齐次方程解公式

$$
\begin{aligned}
u(x,t) &= \frac{1}{2}\left[\varphi(x + at) + \varphi(x - at)\right] + \frac{1}{2a} \int_{x - at}^{x + at} \psi(\xi) d\xi \\
&\quad + \frac{1}{2a} \int_0^t \int_{x - a(t - \tau)}^{x + a(t - \tau)} f(\xi,\tau) d\xi d\tau
\end{aligned} \tag{18}
$$

---

### 3.1.3 半无界问题（延拓法）

#### 情形1：齐次Dirichlet边界

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + f(x,t), & 0 < x < \infty, \; t > 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x < \infty \\
u(0,t) = 0, & t \geq 0
\end{cases} \tag{19}
$$

**解法**：对 $f, \varphi, \psi$ 作**奇延拓**到整个实轴，化为无界问题求解。

#### 情形2：齐次Neumann边界

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + f(x,t), & 0 < x < \infty, \; t > 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x < \infty \\
u_x(0,t) = 0, & t \geq 0
\end{cases} \tag{22}
$$

**解法**：对 $f, \varphi, \psi$ 作**偶延拓**到整个实轴。

#### 情形3：非齐次边界（反射法）

$$
\begin{cases}
u_{tt} = a^2 u_{xx}, & 0 < x < \infty, \; t > 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x < \infty \\
u(0,t) = h(t), & t \geq 0
\end{cases} \tag{24}
$$

**解法**：分区域讨论，考虑特征线 $x = at$。

当 $x \geq at$ 时（波未到达边界）：
$$
u(x,t) = \frac{1}{2}\left[\varphi(x + at) + \varphi(x - at)\right] + \frac{1}{2a} \int_{x - at}^{x + at} \psi(s) ds \tag{25}
$$

当 $0 < x < at$ 时（波已反射）：
$$
u(x,t) = h\left(t - \frac{x}{a}\right) + \frac{1}{2}\left[\varphi(x + at) - \varphi(at - x)\right] + \frac{1}{2a} \int_{at - x}^{x + at} \psi(s) ds \tag{29}
$$

---

## 3.2 三维初值问题（球面平均法、泊松公式）

考虑三维波动方程初值问题：

$$
\begin{cases}
u_{tt} = a^2 (u_{xx} + u_{yy} + u_{zz}), & (x,y,z) \in \mathbb{R}^3, \; t > 0 \\
u|_{t=0} = \varphi(x,y,z), \; u_t|_{t=0} = \psi(x,y,z), & (x,y,z) \in \mathbb{R}^3
\end{cases} \tag{34}
$$

### 3.2.1 球对称解

若解是球对称的：$u(x,y,z;t) = u(\sqrt{x^2 + y^2 + z^2}; t) = u(r;t)$

方程化为：
$$
u_{tt} = a^2 \left(u_{rr} + \frac{2}{r} u_r\right)
$$

两边乘以 $r$ 得：
$$
(ru)_{tt} = a^2 (ru)_{rr}
$$

这是关于 $ru$ 的一维波动方程，其解为：
$$
ru = f_1(r + at) + f_2(r - at)
$$

因此球对称解为：
$$
u(r,t) = \frac{1}{r} \left[f_1(r + at) + f_2(r - at)\right]
$$

---

### 3.2.2 泊松公式（Poisson）

三维波动方程初值问题的解为：

$$
u(M,t) = \frac{1}{4\pi a^2} \frac{\partial}{\partial t} \iint_{S_{at}^M} \frac{\varphi(\xi,\eta,\zeta)}{t} dS + \frac{1}{4\pi a^2} \iint_{S_{at}^M} \frac{\psi(\xi,\eta,\zeta)}{t} dS
$$

其中 $S_{at}^M$ 表示以 $M(x,y,z)$ 为中心，$at$ 为半径的球面。

**球坐标形式**：

$$
\begin{aligned}
u(x,y,z;t) &= \frac{1}{4\pi} \frac{\partial}{\partial t} \int_0^{2\pi} \int_0^\pi \varphi(x + at\sin\theta\cos\phi, \; y + at\sin\theta\sin\phi, \; z + at\cos\theta) \; t\sin\theta \, d\theta d\phi \\
&\quad + \frac{1}{4\pi} \int_0^{2\pi} \int_0^\pi \psi(x + at\sin\theta\cos\phi, \; y + at\sin\theta\sin\phi, \; z + at\cos\theta) \; t\sin\theta \, d\theta d\phi
\end{aligned}
$$

> 💡 **物理意义**：三维空间中波的传播有清晰的波前和波后，即**惠更斯原理**成立，无后效现象。

---

### 3.2.3 三维非齐次方程（推迟势）

$$
\begin{cases}
u_{tt} = a^2 (u_{xx} + u_{yy} + u_{zz}) + f(x,y,z;t), & (x,y,z) \in \mathbb{R}^3, \; t > 0 \\
u|_{t=0} = 0, \; u_t|_{t=0} = 0, & (x,y,z) \in \mathbb{R}^3
\end{cases}
$$

解为推迟势：

$$
u(x,y,z;t) = \frac{1}{4\pi a^2} \iiint_{r \leq at} \frac{f\left(\xi,\eta,\zeta;t - \frac{r}{a}\right)}{r} dV
$$

其中 $r = \sqrt{(x - \xi)^2 + (y - \eta)^2 + (z - \zeta)^2}$。

---

## 3.3 二维初值问题（降维法）

考虑二维波动方程初值问题：

$$
\begin{cases}
u_{tt} = a^2 (u_{xx} + u_{yy}), & (x,y) \in \mathbb{R}^2, \; t > 0 \\
u|_{t=0} = \varphi(x,y), \; u_t|_{t=0} = \psi(x,y), & (x,y) \in \mathbb{R}^2
\end{cases} \tag{37}
$$

### 降维法

将二维问题视为三维问题，其初始数据与 $z$ 无关。代入三维泊松公式化简得：

$$
\begin{aligned}
u(x,y;t) &= \frac{1}{2\pi a} \frac{\partial}{\partial t} \iint_{\Sigma_{at}^M} \frac{\varphi(\xi,\eta) d\xi d\eta}{\sqrt{a^2 t^2 - (\xi - x)^2 - (\eta - y)^2}} \\
&\quad + \frac{1}{2\pi a} \iint_{\Sigma_{at}^M} \frac{\psi(\xi,\eta) d\xi d\eta}{\sqrt{a^2 t^2 - (\xi - x)^2 - (\eta - y)^2}}
\end{aligned}
$$

其中 $\Sigma_{at}^M$ 表示以 $M(x,y)$ 为中心，$at$ 为半径的圆盘。

**极坐标形式**：

$$
\begin{aligned}
u(x,y;t) &= \frac{1}{2\pi a} \frac{\partial}{\partial t} \int_0^{at} \int_0^{2\pi} \frac{\varphi(x + r\cos\theta, y + r\sin\theta) r d\theta dr}{\sqrt{a^2 t^2 - r^2}} \\
&\quad + \frac{1}{2\pi a} \int_0^{at} \int_0^{2\pi} \frac{\psi(x + r\cos\theta, y + r\sin\theta) r d\theta dr}{\sqrt{a^2 t^2 - r^2}}
\end{aligned}
$$

> 💡 **物理意义**：二维空间中波的传播有**弥散现象**，即波过后有后效，不像三维有清晰的波前和波后。

---

## 3.4 依赖区域、决定区域、影响区域

### 依赖区域

点 $(x,t)$ 处的解 $u(x,t)$ 只依赖于初始轴上区间 $[x - at, x + at]$ 内的初值，这个区间称为点 $(x,t)$ 的**依赖区域**。

### 决定区域

初始轴上区间 $[x_1, x_2]$ 的初值唯一决定了区域：
$$
D = \{(x,t) \mid x_1 + at \leq x \leq x_2 - at, \; t \geq 0\}
$$
内各点的解值，$D$ 称为区间 $[x_1, x_2]$ 的**决定区域**。

### 影响区域

初始轴上点 $x_0$ 处的初扰动能影响到区域：
$$
G = \{(x,t) \mid x_0 - at \leq x \leq x_0 + at, \; t \geq 0\}
$$
内各点的解值，$G$ 称为点 $x_0$ 的**影响区域**。

> 💡 **特征线的意义**：特征线 $x \pm at = C$ 是波的传播路径，也是扰动的分界线。

---

## 📝 考试真题与解题模板

### 达朗贝尔公式解题步骤（必背）

| 步骤 | 操作 |
|------|------|
| **Step 1** | 识别问题类型：齐次/非齐次，无界/半无界 |
| **Step 2** | 无界齐次问题直接套用达朗贝尔公式 |
| **Step 3** | 非齐次问题用齐次化原理叠加 |
| **Step 4** | 半无界问题用延拓法或反射法 |
| **Step 5** | 计算积分，注意上下限 |

---

### 真题1（2023级期末试题）⭐⭐⭐⭐⭐

**题目**：利用达朗贝尔公式求解波动方程初值问题

$$
\begin{cases}
u_{tt} = u_{xx} + t\cos x, & -\infty < x < +\infty, \; t > 0 \\
u(x,0) = 1 + x^2, \quad u_t(x,0) = \cos x, & -\infty < x < +\infty
\end{cases}
$$

**解答：**

这是非齐次方程初值问题，用叠加原理。

**Step 1. 齐次部分解 $v$**

$\varphi(x) = 1 + x^2, \quad \psi(x) = \cos x$

由达朗贝尔公式：
$$
v(x,t) = \frac{1}{2}\left[\varphi(x + t) + \varphi(x - t)\right] + \frac{1}{2} \int_{x - t}^{x + t} \cos s \, ds
$$

计算：
$$
\varphi(x + t) + \varphi(x - t) = [1 + (x + t)^2] + [1 + (x - t)^2] = 2 + 2x^2 + 2t^2
$$

$$
\int_{x - t}^{x + t} \cos s \, ds = \sin(x + t) - \sin(x - t) = 2\cos x \sin t
$$

因此：
$$
v(x,t) = 1 + x^2 + t^2 + \cos x \sin t
$$

**Step 2. 非齐次部分解 $w$（齐次化原理）**

$$
w(x,t) = \frac{1}{2} \int_0^t \int_{x - (t - \tau)}^{x + (t - \tau)} \tau \cos \xi \, d\xi d\tau
$$

内层积分：
$$
\int_{x - (t - \tau)}^{x + (t - \tau)} \cos \xi \, d\xi = 2\cos x \sin(t - \tau)
$$

外层积分：
$$
w(x,t) = \int_0^t \tau \cos x \sin(t - \tau) d\tau = \cos x \int_0^t \tau \sin(t - \tau) d\tau
$$

分部积分：
$$
\int_0^t \tau \sin(t - \tau) d\tau = t - \sin t
$$

因此：
$$
w(x,t) = \cos x (t - \sin t)
$$

**Step 3. 最终解**
$$
u(x,t) = v(x,t) + w(x,t) = 1 + x^2 + t^2 + t\cos x
$$

---

### 真题2（2020级期末试题）⭐⭐⭐⭐⭐

**题目**：一维波动方程初值问题

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + \cos(x + at), & -\infty < x < +\infty, \; t > 0 \\
u(x,0) = x, \quad u_t(x,0) = \sin x, & -\infty < x < +\infty
\end{cases}
$$

**解答：**

**Step 1. 齐次部分解 $v$**

$$
v(x,t) = \frac{1}{2}\left[(x + at) + (x - at)\right] + \frac{1}{2a} \int_{x - at}^{x + at} \sin s \, ds
$$

$$
= x + \frac{1}{2a} \left[-\cos(x + at) + \cos(x - at)\right]
$$

$$
= x + \frac{1}{a} \sin x \sin at
$$

**Step 2. 非齐次部分解 $w$**

用齐次化原理：
$$
w(x,t) = \frac{1}{2a} \int_0^t \int_{x - a(t - \tau)}^{x + a(t - \tau)} \cos(\xi + a\tau) d\xi d\tau
$$

内层积分：
$$
\int_{x - a(t - \tau)}^{x + a(t - \tau)} \cos(\xi + a\tau) d\xi = \sin(x + at) - \sin(x - at + 2a\tau)
$$

外层积分：
$$
w(x,t) = \frac{1}{2a} \int_0^t \left[\sin(x + at) - \sin(x - at + 2a\tau)\right] d\tau
$$

$$
= \frac{t}{2a} \sin(x + at) + \frac{1}{4a^2} \left[\cos(x + at) - \cos(x - at)\right]
$$

**Step 3. 最终解**

$$
\begin{aligned}
u(x,t) &= x + \frac{1}{a} \sin x \sin at + \frac{t}{2a} \sin(x + at) \\
&\quad + \frac{1}{4a^2} \left[\cos(x + at) - \cos(x - at)\right]
\end{aligned}
$$

---

## 本章小结

| 方法 | 适用问题 | 核心公式 | 考试频率 |
|------|----------|----------|----------|
| **达朗贝尔公式** | 一维齐次无界问题 | $u = \frac{1}{2}[\varphi(x+at)+\varphi(x-at)] + \frac{1}{2a}\int_{x-at}^{x+at}\psi(s)ds$ | ⭐⭐⭐⭐⭐ |
| **齐次化原理** | 一维非齐次问题 | $w = \frac{1}{2a}\int_0^t\int_{x-a(t-\tau)}^{x+a(t-\tau)}f(\xi,\tau)d\xi d\tau$ | ⭐⭐⭐⭐ |
| **延拓法** | 半无界问题 | 奇延拓（$u=0$）、偶延拓（$u_x=0$） | ⭐⭐⭐ |
| **泊松公式** | 三维问题 | 球面积分 | ⭐⭐⭐ |
| **降维法** | 二维问题 | 圆域积分 | ⭐⭐ |

> 💡 **考试提醒**：达朗贝尔公式每年必考！重点掌握齐次化原理的应用和非齐次项的积分计算。

---

> 本文件基于PPT内容整理优化，补充了历年考试真题和详细解题过程。

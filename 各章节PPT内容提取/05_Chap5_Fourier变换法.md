# 第5章 Fourier变换法

> ⭐ **考试重点**：Fourier变换性质、热传导方程、Laplace方程为高频考点，常考解答题

---

## 目录

- [5.1 Fourier积分定理](#51-fourier积分定理)
- [5.2 Fourier变换及其性质](#52-fourier变换及其性质)
- [5.3 Fourier变换应用](#53-fourier变换应用)
- [5.4 正弦/余弦变换](#54-正弦余弦变换)
- [📝 考试真题与解题模板](#-考试真题与解题模板)

---

## 5.1 Fourier积分定理

### Fourier积分定理

**定理**：若函数 $f(x)$ 在 $(-\infty, +\infty)$ 上分段光滑且绝对可积，则在 $f(x)$ 的连续点处，成立：

$$
f(x) = \frac{1}{\pi} \int_{0}^{+\infty} \int_{-\infty}^{+\infty} f(t) \cos\lambda(t - x) dt d\lambda \tag{1}
$$

在间断点处，上式右端收敛于 $\frac{1}{2}(f(x+0) + f(x-0))$。

---

### Fourier变换定义

定义 $f(x)$ 的 **Fourier变换** 为：

$$
F(\lambda) = \mathcal{F}[f(x)] = \int_{-\infty}^{+\infty} f(t) e^{-i\lambda t} dt \tag{2}
$$

则 **Fourier逆变换** 为：

$$
f(x) = \mathcal{F}^{-1}[F(\lambda)] = \frac{1}{2\pi} \int_{-\infty}^{+\infty} F(\lambda) e^{i\lambda x} d\lambda \tag{3}
$$

其中 $F(\lambda)$ 称为 $f(x)$ 的 **象函数**，$f(x)$ 称为 $F(\lambda)$ 的 **象原函数**。

---

## 5.2 Fourier变换及其性质

### 正弦变换与余弦变换

#### 正弦变换（奇函数延拓）

若 $f(x)$ 是奇函数，则有 Fourier 正弦积分公式：

$$
f(x) = \frac{2}{\pi} \int_{0}^{+\infty} \sin\lambda x \left( \int_{0}^{+\infty} f(t) \sin\lambda t dt \right) d\lambda \tag{4}
$$

定义 **Fourier正弦变换**：

$$
F_s(\lambda) = \int_{0}^{+\infty} f(t) \sin\lambda t dt \tag{5}
$$

则正弦逆变换为：

$$
f(x) = \frac{2}{\pi} \int_{0}^{+\infty} F_s(\lambda) \sin\lambda x d\lambda \tag{6}
$$

> 💡 **适用场景**：具有 Dirichlet 边界条件 $u(0,t) = 0$ 的半无界问题

---

#### 余弦变换（偶函数延拓）

若 $f(x)$ 是偶函数，则有 Fourier 余弦积分公式：

$$
f(x) = \frac{2}{\pi} \int_{0}^{+\infty} \cos\lambda x \left( \int_{0}^{+\infty} f(t) \cos\lambda t dt \right) d\lambda \tag{7}
$$

定义 **Fourier余弦变换**：

$$
F_c(\lambda) = \int_{0}^{+\infty} f(t) \cos\lambda t dt \tag{8}
$$

则余弦逆变换为：

$$
f(x) = \frac{2}{\pi} \int_{0}^{+\infty} F_c(\lambda) \cos\lambda x d\lambda \tag{8}
$$

> 💡 **适用场景**：具有 Neumann 边界条件 $u_x(0,t) = 0$ 的半无界问题

---

### 基本性质

#### 1. 线性性质

Fourier变换是线性变换：

$$
\mathcal{F}[c_1 f + c_2 g] = c_1 \mathcal{F}[f] + c_2 \mathcal{F}[g] = c_1 F(\lambda) + c_2 G(\lambda)
$$

$$
\mathcal{F}^{-1}[c_1 F(\lambda) + c_2 G(\lambda)] = c_1 \mathcal{F}^{-1}[F(\lambda)] + c_2 \mathcal{F}^{-1}[G(\lambda)]
$$

其中 $c_1, c_2$ 为任意常数。

---

#### 2. 相似性质（伸缩性质）

若 $a \neq 0$ 为常数，则：

$$
\mathcal{F}[f(at)] = \frac{1}{|a|} F\left(\frac{\lambda}{a}\right)
$$

---

#### 3. 微分性质

若 $f(t), f'(t), \cdots, f^{(n)}(t)$ 在 $(-\infty, +\infty)$ 上分段连续可微，且当 $|t| \to \infty$ 时，各阶导数趋于0，则：

$$
\mathcal{F}[f'(t)] = i\lambda \mathcal{F}[f(t)]
$$

一般地，对 $n$ 阶导数有：

$$
\mathcal{F}[f^{(n)}(t)] = (i\lambda)^n \mathcal{F}[f(t)]
$$

> 💡 **微分方程降维**：该性质可将偏微分方程转化为常微分方程

---

#### 4. 乘多项式性质

若 $f(t)$ 在 $(-\infty, +\infty)$ 上绝对可积，则：

$$
\mathcal{F}[t^n f(t)] = i^n F^{(n)}(\lambda)
$$

其中 $F^{(n)}(\lambda)$ 是 $F(\lambda)$ 对 $\lambda$ 的 $n$ 阶导数。

---

### 重要例题：高斯函数的Fourier变换

求函数 $f(x) = e^{-bx^2}$ 的 Fourier 变换，其中常数 $b > 0$。

**解**：

$$
F(\lambda) = \int_{-\infty}^{+\infty} e^{-bt^2} e^{-i\lambda t} dt = \sqrt{\frac{\pi}{b}} e^{-\lambda^2/(4b)}
$$

> 💡 **重要结论**：高斯函数的Fourier变换仍是高斯函数！这是一个非常重要的结论，常用于热传导方程求解。

---

### 卷积定理

#### 卷积定义

函数 $f(x)$ 与 $g(x)$ 的卷积定义为：

$$
(f * g)(x) = \int_{-\infty}^{+\infty} f(x - t) g(t) dt \tag{9}
$$

#### 卷积性质

1. **交换律**：$(f * g)(x) = (g * f)(x)$
2. **结合律**：$(f * (g * h))(x) = ((f * g) * h)(x)$
3. **分配律**：$(f * (g + h))(x) = (f * g)(x) + (f * h)(x)$

#### 卷积定理

若 $F(\lambda) = \mathcal{F}[f(x)], G(\lambda) = \mathcal{F}[g(x)]$，则：

$$
\mathcal{F}[(f * g)(x)] = F(\lambda) \cdot G(\lambda)
$$

$$
\mathcal{F}^{-1}[F(\lambda) G(\lambda)] = \frac{1}{2\pi} (f * g)(x)
$$

> 💡 **乘积定理**：象函数的乘积对应象原函数卷积的 $1/(2\pi)$ 倍

---

## 5.3 Fourier变换应用

### 例1：无限长杆热传导方程初值问题

$$
\begin{cases}
u_t = a^2 u_{xx} + f(x, t), & x \in \mathbb{R}, t > 0 \\
u(x, 0) = \varphi(x), & x \in \mathbb{R}^1
\end{cases} \tag{11}
$$

**求解过程**：

设 $U(\lambda, t), F(\lambda, t), \Phi(\lambda)$ 分别为 $u(x,t), f(x,t), \varphi(x)$ 关于 $x$ 的 Fourier 变换。

对方程两边关于 $x$ 作 Fourier 变换，得到常微分方程初值问题：

$$
\begin{cases}
\frac{dU}{dt} + a^2 \lambda^2 U = F(\lambda, t), & t > 0 \\
U(\lambda, 0) = \Phi(\lambda)
\end{cases} \tag{12}
$$

其解为：

$$
U(\lambda, t) = \Phi(\lambda) e^{-a^2 \lambda^2 t} + \int_{0}^{t} F(\lambda, \tau) e^{-a^2 \lambda^2 (t - \tau)} d\tau \tag{13}
$$

作 Fourier 逆变换，利用卷积定理：

$$
\mathcal{F}^{-1}\left[e^{-a^2 \lambda^2 t}\right] = \frac{1}{2a\sqrt{\pi t}} \exp\left(-\frac{x^2}{4a^2 t}\right)
$$

这就是热传导方程的基本解——**高斯核函数**。

由卷积定理：

$$
\mathcal{F}^{-1}\left[\Phi(\lambda) e^{-a^2 \lambda^2 t}\right] = \mathcal{F}^{-1}[\Phi(\lambda)] * \mathcal{F}^{-1}\left[e^{-a^2 \lambda^2 t}\right]
$$

因此，原问题的解为：

$$
\begin{aligned}
u(x, t) &= \frac{1}{2a\sqrt{\pi t}} \int_{-\infty}^{+\infty} \varphi(\xi) \exp\left(-\frac{(x - \xi)^2}{4a^2 t}\right) d\xi \\
&\quad + \int_{0}^{t} \int_{-\infty}^{+\infty} \frac{f(\xi, \tau)}{2a\sqrt{\pi (t - \tau)}} \exp\left(-\frac{(x - \xi)^2}{4a^2 (t - \tau)}\right) d\xi d\tau
\end{aligned} \tag{15}
$$

---

### 例2：上半平面Laplace方程Dirichlet问题

考虑上半平面 $y > 0$ 的调和函数：

$$
\begin{cases}
u_{xx} + u_{yy} = 0, & x \in \mathbb{R}^1, y > 0 \\
u(x, 0) = f(x), & x \in \mathbb{R}^1 \\
\lim\limits_{|x| \to \infty} u(x, y) = 0, \lim\limits_{|x| \to \infty} u_x(x, y) = 0
\end{cases} \tag{16}
$$

**求解过程**：

设 $U(\lambda, y), F(\lambda)$ 分别为 $u(x, y), f(x)$ 关于 $x$ 的 Fourier 变换。

对方程两边作 Fourier 变换：

$$
\frac{d^2 U}{dy^2} - \lambda^2 U = 0, \quad U(\lambda, 0) = F(\lambda) \tag{17}
$$

其通解为：

$$
U(\lambda, y) = A(\lambda) e^{-\lambda y} + B(\lambda) e^{\lambda y} \tag{18}
$$

由有界性条件 $y \to \infty$ 时 $u$ 有界，得：
- 当 $\lambda > 0$ 时，$B(\lambda) = 0$
- 当 $\lambda < 0$ 时，$A(\lambda) = 0$

因此：

$$
U(\lambda, y) = F(\lambda) e^{-|\lambda| y}
$$

作逆变换：

$$
u(x, y) = \frac{1}{2\pi} \int_{-\infty}^{+\infty} F(\lambda) e^{-|\lambda| y} e^{i\lambda x} d\lambda \tag{19}
$$

利用卷积定理，可得：

$$
u(x, y) = \frac{y}{\pi} \int_{-\infty}^{+\infty} \frac{f(\xi)}{(x - \xi)^2 + y^2} d\xi \tag{20}
$$

这就是著名的 **Poisson积分公式**。

---

### 例3：无限长弦振动方程初值问题

$$
\begin{cases}
u_{tt} = a^2 u_{xx}, & x \in \mathbb{R}^1, t > 0 \\
u(x, 0) = f(x), u_t(x, 0) = g(x), & x \in \mathbb{R}^1
\end{cases} \tag{21}
$$

**求解过程**：

设 $U(\lambda, t), F(\lambda), G(\lambda)$ 分别为 $u(x, t), f(x), g(x)$ 关于 $x$ 的 Fourier 变换。

对方程两边作 Fourier 变换：

$$
\begin{cases}
\frac{d^2 U}{dt^2} + a^2 \lambda^2 U = 0, & t > 0 \\
U(\lambda, 0) = F(\lambda), U_t(\lambda, 0) = G(\lambda)
\end{cases} \tag{22}
$$

其解为：

$$
U(\lambda, t) = F(\lambda) \cos a\lambda t + \frac{G(\lambda)}{a\lambda} \sin a\lambda t
$$

作逆变换，利用三角恒等式化简：

$$
\begin{aligned}
u(x, t) &= \frac{1}{4\pi} \int_{-\infty}^{+\infty} F(\lambda) \left(e^{i\lambda (x + at)} + e^{i\lambda (x - at)}\right) d\lambda \\
&\quad + \frac{1}{4\pi a} \int_{-\infty}^{+\infty} \frac{G(\lambda)}{i\lambda} \left(e^{i\lambda (x + at)} - e^{i\lambda (x - at)}\right) d\lambda
\end{aligned}
$$

最终得到达朗贝尔公式：

$$
u(x, t) = \frac{1}{2}\left[f(x + at) + f(x - at)\right] + \frac{1}{2a} \int_{x - at}^{x + at} g(\xi) d\xi \tag{24}
$$

---

## 5.4 正弦/余弦变换

### 半无界问题的应用

#### 例4：用Fourier余弦变换求解

$$
\begin{cases}
u_t = a^2 u_{xx}, & x > 0, t > 0 \\
u(x, 0) = 0, u(0, t) = u_0, & t \geq 0 \\
\lim\limits_{x \to \infty} u(x, t) = 0, \lim\limits_{x \to \infty} u_x(x, t) = 0
\end{cases} \tag{25}
$$

#### 例5：用Fourier余弦变换求解

$$
\begin{cases}
u_t = a^2 u_{xx}, & x > 0, t > 0 \\
u(x, 0) = 0, u_x(0, t) = u_0, & t \geq 0 \\
\lim\limits_{x \to \infty} u(x, t) = 0, \lim\limits_{x \to \infty} u_x(x, t) = 0
\end{cases} \tag{26}
$$

**解法提示**：由于边界条件是 $u_x(0,t)$ 已知，适合使用Fourier余弦变换。

设 $U_c(\lambda, t) = \mathcal{F}_c[u(x, t)] = \int_{0}^{+\infty} u(x, t) \cos\lambda x dx$

对方程两边关于 $x$ 作Fourier余弦变换：

$$
\frac{dU_c}{dt} + a^2 \lambda^2 U_c = a^2 u_0, \quad t > 0
$$

其解为：

$$
U_c(\lambda, t) = \frac{u_0}{\lambda^2} \left(1 - e^{-a^2 \lambda^2 t}\right) \tag{27}
$$

作Fourier余弦逆变换：

$$
u(x, t) = \mathcal{F}_c^{-1}[U_c(\lambda, t)] = \frac{2u_0}{\pi} \int_{0}^{+\infty} \frac{\cos\lambda x}{\lambda^2} \left(1 - e^{-a^2 \lambda^2 t}\right) d\lambda
$$

---

## 📝 考试真题与解题模板

### Fourier变换法解题步骤（必背）

| 步骤 | 操作 |
|------|------|
| **Step 1** | 对定解问题关于空间变量作Fourier变换（选择合适类型：正/弦/余弦） |
| **Step 2** | 原偏微分方程转化为关于时间的常微分方程 |
| **Step 3** | 对初始条件作相应的Fourier变换 |
| **Step 4** | 求解常微分方程初值问题，得到象函数 $U(\lambda, t)$ |
| **Step 5** | 作Fourier逆变换，利用卷积定理得到原解 |

---

### 真题1：热传导方程初值问题

**题目**：利用Fourier变换方法，求解下列热传导方程初值问题

$$
\begin{cases}
u_t = a^2 u_{xx} + t^2, & -\infty < x < +\infty, \; t > 0 \\
u(x, 0) = x^2 + \sin x, & -\infty < x < +\infty
\end{cases}
$$

已知 $\int_{0}^{+\infty} e^{-a^2 x^2} \cos bx dx = \frac{\sqrt{\pi}}{2a} e^{-b^2/(4a^2)}$

---

## 本章小结

| 类型 | 公式表示 | 考试频率 |
|------|----------|----------|
| **Fourier变换** | $F(\lambda) = \int_{-\infty}^{+\infty} f(t) e^{-i\lambda t} dt$ | ⭐⭐⭐⭐⭐ |
| **微分性质** | $\mathcal{F}[f^{(n)}(t)] = (i\lambda)^n F(\lambda)$ | ⭐⭐⭐⭐⭐ |
| **卷积定理** | $\mathcal{F}[f * g] = F(\lambda) G(\lambda)$ | ⭐⭐⭐⭐ |
| **热传导方程解** | $u = \varphi * \frac{1}{2a\sqrt{\pi t}} e^{-x^2/(4a^2 t)}$ | ⭐⭐⭐⭐⭐ |
| **Poisson公式** | $u(x,y) = \frac{y}{\pi} \int_{-\infty}^{+\infty} \frac{f(\xi)}{(x-\xi)^2+y^2} d\xi$ | ⭐⭐⭐⭐ |

> 💡 **考试提醒**：Fourier变换法每年必考！重点掌握热传导方程求解、卷积定理应用，注意区分无界问题和半无界问题（正变换 vs 正弦/余弦变换）。

---

> 本文件基于PPT内容整理优化，补充了详细的公式推导和解题步骤。

# 第6章 Laplace变换法

> ⭐ **考试重点**：微分性质、卷积定理、常微分方程求解为高频考点，偏微分方程应用为解答题常考内容

---

## 目录

- [6.1 Laplace变换定义与存在性](#61-laplace变换定义与存在性)
- [6.2 Laplace变换基本性质](#62-laplace变换基本性质)
- [6.3 卷积定理](#63-卷积定理)
- [6.4 Laplace逆变换](#64-laplace逆变换)
- [6.5 应用举例](#65-应用举例)
- [📝 考试真题与解题模板](#-考试真题与解题模板)

---

## 6.1 Laplace变换定义与存在性

### Laplace变换定义

**定义**：设函数 $f(t)$ 定义在 $[0, +\infty)$ 上，若含参变量的无穷积分

$$
\int_{0}^{+\infty} f(t) e^{-st} dt
$$

收敛，则称此积分为函数 $f(t)$ 的 **Laplace变换**，记为：

$$
\mathcal{L}[f(t)] = F(s) = \int_{0}^{+\infty} f(t) e^{-st} dt
$$

相应地，**Laplace逆变换** 定义为：

$$
f(t) = \mathcal{L}^{-1}[F(s)] = \frac{1}{2\pi i} \int_{\sigma - i\infty}^{\sigma + i\infty} F(s) e^{st} ds \tag{1}
$$

其中常数 $\sigma > 0$ 充分大，上式为复围道积分，可用留数定理计算。

---

### Laplace变换存在定理

**定理**：若函数 $f(t)$ 定义在 $[0, +\infty)$ 上，满足下列条件：

1. **(L1)** 在 $t \geq 0$ 的任一有限区间上分段连续；
2. **(L2)** 当 $t \to +\infty$ 时，$f(t)$ 的增长是指数级的，即存在常数 $M > 0, a \geq 0$，使得

$$
|f(t)| \leq M e^{at}, \quad t \geq 0 \tag{2}
$$

则称 $f(t)$ 是指数级增长的，$a$ 称为增长指数。

此时，$f(t)$ 的 Laplace 变换 $F(s)$ 在半平面 $\text{Re}(s) > a$ 上存在且解析。

> 💡 **记忆**：满足上述条件的函数记为 $f(t) = O(e^{at})$。

---

## 6.2 Laplace变换基本性质

### 1. 线性性质

Laplace变换是线性变换：

$$
\mathcal{L}[c_1 f(t) + c_2 g(t)] = c_1 \mathcal{L}[f(t)] + c_2 \mathcal{L}[g(t)] = c_1 F(s) + c_2 G(s)
$$

$$
\mathcal{L}^{-1}[c_1 F(s) + c_2 G(s)] = c_1 \mathcal{L}^{-1}[F(s)] + c_2 \mathcal{L}^{-1}[G(s)]
$$

其中 $c_1, c_2$ 为任意常数。

---

### 2. 微分性质（最重要！）

若 $f'(t)$ 在 $[0, +\infty)$ 上分段连续，$f(t)$ 在 $[0, +\infty)$ 上连续，且当 $t \to \infty$ 时是指数级的，则

$$
\mathcal{L}[f'(t)] = s F(s) - f(0)
$$

**推广到n阶导数**：

若 $f(t), f'(t), \cdots, f^{(n-1)}(t)$ 在 $[0, +\infty)$ 上连续，$f^{(n)}(t)$ 在 $[0, +\infty)$ 上分段连续，则

$$
\begin{aligned}
\mathcal{L}[f^{(n)}(t)] &= s^n F(s) - s^{n-1} f(0) - s^{n-2} f'(0) - \cdots \\
&\quad - s f^{(n-2)}(0) - f^{(n-1)}(0)
\end{aligned}
$$

> 💡 **微分方程求解关键**：该性质可将常微分方程转化为代数方程！

---

### 3. 积分性质

$$
\mathcal{L}^{-1}\left[\frac{F(s)}{s}\right] = \int_{0}^{t} f(\tau) d\tau
$$

**推广**：

$$
\mathcal{L}^{-1}\left[\frac{F(s)}{s^n}\right] = \int_{0}^{t} \int_{0}^{t_{n-1}} \cdots \int_{0}^{t_1} f(\tau) d\tau dt_1 \cdots dt_{n-1}
$$

对 $n = 1, 2, \cdots$ 成立。

---

### 4. 乘多项式性质（象函数微分）

$$
\mathcal{L}[t^n f(t)] = (-1)^n \frac{d^n}{ds^n} F(s)
$$

---

### 5. 位移性质

第一类位移性质（象函数位移）：

$$
\mathcal{L}[e^{at} f(t)] = F(s - a)
$$

第二类位移性质（延迟性质）：

$$
\mathcal{L}^{-1}\left[e^{-as} F(s)\right] = u_a(t) f(t - a)
$$

其中 $a \geq 0$，且
$$
u_a(t) = H(t - a) =
\begin{cases}
0, & t < a \\
1, & t \geq a
\end{cases} \tag{3}
$$

称为Heaviside单位阶跃函数。

---

### 6. 象函数积分性质

若 $\mathcal{L}[f(t)] = F(s)$，且积分 $\int_{s}^{\infty} F(\sigma) d\sigma$ 收敛，则

$$
\mathcal{L}\left[\frac{f(t)}{t}\right] = \int_{s}^{\infty} F(\sigma) d\sigma
$$

**推导**：交换积分顺序

$$
\begin{aligned}
\int_{s}^{\infty} F(\sigma) d\sigma &= \int_{s}^{\infty} \int_{0}^{\infty} e^{-\sigma t} f(t) dt d\sigma \\
&= \int_{0}^{\infty} \int_{s}^{\infty} f(t) e^{-\sigma t} d\sigma dt \\
&= \int_{0}^{\infty} \frac{f(t)}{t} e^{-st} dt = \mathcal{L}\left[\frac{f(t)}{t}\right]
\end{aligned}
$$

---

## 6.3 卷积定理

### 卷积定义

函数 $f(t)$ 与 $g(t)$ 在 $[0, +\infty)$ 上的卷积定义为：

$$
(f * g)(t) = \int_{0}^{t} f(t - \tau) g(\tau) d\tau
$$

> 💡 **注意**：与Fourier卷积的区别是积分限不同！

---

### 卷积定理

若 $\mathcal{L}[f(t)] = F(s), \mathcal{L}[g(t)] = G(s)$，则

$$
\mathcal{L}[(f * g)(t)] = \mathcal{L}[f(t)] \mathcal{L}[g(t)] = F(s) G(s)
$$

或

$$
\mathcal{L}^{-1}[F(s) G(s)] = (f * g)(t)
$$

> 💡 **核心应用**：乘积的逆变换等于象原函数的卷积，这是求解非齐次项的关键工具！

---

## 6.4 Laplace逆变换

### 基本变换表（必须熟记！）

| 原函数 $f(t)$ | 象函数 $F(s)$ | 备注 |
|---------------|---------------|------|
| $1$ | $\frac{1}{s}$ | $s > 0$ |
| $t^n$ | $\frac{n!}{s^{n+1}}$ | $n = 0, 1, 2, \cdots$ |
| $e^{at}$ | $\frac{1}{s - a}$ | $s > a$ |
| $\cos bt$ | $\frac{s}{s^2 + b^2}$ | $s > 0$ |
| $\sin bt$ | $\frac{b}{s^2 + b^2}$ | $s > 0$ |
| $e^{at} \cos bt$ | $\frac{s - a}{(s - a)^2 + b^2}$ | |
| $e^{at} \sin bt$ | $\frac{b}{(s - a)^2 + b^2}$ | |

---

### 常用计算技巧

#### 例1：利用乘多项式性质

$$
\mathcal{L}[t \sin 3t] = -\frac{d}{ds} \left(\frac{3}{s^2 + 9}\right) = \frac{6s}{(s^2 + 9)^2}
$$

利用位移性质：

$$
\mathcal{L}[e^{-4t} t \sin 3t] = \frac{6(s + 4)}{((s + 4)^2 + 9)^2}
$$

---

#### 例2：利用三角恒等式化简

$$
\sin^3 2t = \frac{1}{2}(1 - \cos 4t) \sin 2t = \frac{1}{2}(\sin 2t - \cos 4t \sin 2t)
$$

$$
= \frac{3}{4} \sin 2t - \frac{1}{4} \sin 6t
$$

因此

$$
\begin{aligned}
\mathcal{L}[\sin^3 2t] &= \frac{3}{4} \mathcal{L}[\sin 2t] - \frac{1}{4} \mathcal{L}[\sin 6t] \\
&= \frac{3}{4} \cdot \frac{2}{s^2 + 4} - \frac{1}{4} \cdot \frac{6}{s^2 + 36} \\
&= \frac{48}{(s^2 + 4)(s^2 + 36)}
\end{aligned}
$$

---

#### 例3：利用象函数积分性质

$$
\mathcal{L}[1 - \cos t] = \frac{1}{s} - \frac{s}{s^2 + 1}
$$

因此

$$
\begin{aligned}
\mathcal{L}\left[\frac{1 - \cos t}{t}\right] &= \int_{s}^{\infty} \left(\frac{1}{\sigma} - \frac{\sigma}{\sigma^2 + 1}\right) d\sigma \\
&= \left.\ln \frac{\sigma}{\sqrt{\sigma^2 + 1}}\right|_{s}^{\infty} = \frac{1}{2} \ln \frac{s^2 + 1}{s^2}
\end{aligned}
$$

---

### 逆变换举例

#### 例1：利用延迟性质

求 $F(s) = \frac{e^{-as}}{1 + s^2}, a > 0$ 的逆变换。

解：由于

$$
\mathcal{L}^{-1}\left[\frac{1}{1 + s^2}\right] = \sin t
$$

由延迟性质得

$$
f(t) = \mathcal{L}^{-1}[F(s)] = \sin(t - a) H(t - a)
$$

---

#### 例2：部分分式分解

求 $F(s) = \frac{1}{(s^2 + \omega^2)^2}$ 的逆变换。

解：利用代数变形

$$
\frac{1}{(s^2 + \omega^2)^2} = \frac{1}{2\omega^2} \left( \frac{1}{s^2 + \omega^2} - \frac{s^2 - \omega^2}{(s^2 + \omega^2)^2} \right)
$$

已知

$$
\mathcal{L}^{-1}\left[\frac{1}{s^2 + \omega^2}\right] = \frac{\sin \omega t}{\omega}
$$

$$
\mathcal{L}^{-1}\left[\frac{s^2 - \omega^2}{(s^2 + \omega^2)^2}\right] = t \cos \omega t
$$

因此

$$
f(t) = \frac{1}{2\omega^3} [\sin \omega t - \omega t \cos \omega t]
$$

---

## 6.5 应用举例

### Laplace变换法解题步骤（必背）

| 步骤 | 操作 |
|------|------|
| **Step 1** | 对定解问题关于时间变量 $t$ 作Laplace变换 |
| **Step 2** | 原偏微分方程转化为关于空间变量的常微分方程 |
| **Step 3** | 对边界条件作相应的Laplace变换 |
| **Step 4** | 求解常微分方程边值问题，得到象函数 $U(x, s)$ |
| **Step 5** | 作Laplace逆变换，得到原解 |

---

### 例1：常微分方程初值问题

求解

$$
y''(t) + 4y'(t) + 3y(t) = 0, \quad t > 0
$$

初始条件：$y(0) = 3, y'(0) = 1$

**解**：设 $Y(s) = \mathcal{L}[y(t)]$，对方程两边作Laplace变换

$$
(s^2 Y(s) - s y(0) - y'(0)) + 4(s Y(s) - y(0)) + 3Y(s) = 0
$$

整理得

$$
(s + 3)(s + 1)Y(s) = 3s + 13
$$

$$
Y(s) = \frac{3s + 13}{(s + 1)(s + 3)} = \frac{5}{s + 1} - \frac{2}{s + 3}
$$

作逆变换

$$
y(t) = \mathcal{L}^{-1}[Y(s)] = 5e^{-t} - 2e^{-3t}
$$

---

### 例2：半无界弦振动方程

$$
\begin{cases}
u_{tt} = c^2 u_{xx} + f(t), & 0 < x < \infty, t > 0 \\
u(x, 0) = 0, u_t(x, 0) = 0, & 0 \leq x < \infty \\
u(0, t) = 0, \lim\limits_{x \to \infty} u_x(x, t) = 0
\end{cases} \tag{4}
$$

**解**：设 $U(x, s) = \mathcal{L}[u(x, t)], F(s) = \mathcal{L}[f(t)]$，则变换后方程为

$$
\frac{d^2 U}{dx^2} - \frac{s^2}{c^2} U = -\frac{F(s)}{c^2} \tag{5}
$$

边界条件变换后：$U(0, s) = 0, \lim\limits_{x \to \infty} U_x(x, s) = 0$

通解为

$$
U(x, s) = A e^{sx/c} + B e^{-sx/c}
$$

由边界条件得 $A = 0, B = \frac{F(s)}{s^2}$，因此

$$
U(x, s) = \frac{F(s)}{s^2} (1 - e^{-sx/c})
$$

作逆变换，利用卷积定理

$$
\mathcal{L}^{-1}\left[\frac{F(s)}{s^2}\right] = \int_{0}^{t} \int_{0}^{\tau} f(\sigma) d\sigma d\tau = G(t)
$$

由延迟性质

$$
\mathcal{L}^{-1}\left[\frac{F(s) e^{-sx/c}}{s^2}\right] = H(t - x/c) G(t - x/c)
$$

因此

$$
u(x, t) = G(t) - H(t - x/c) G(t - x/c)
$$

其中 $H(t - a)$ 是Heaviside单位阶跃函数。

---

### 例3：有界弦振动方程

$$
\begin{cases}
u_{tt} = u_{xx}, & 0 < x < 1, t > 0 \\
u(0, t) = u(1, t) = 0, & t \geq 0 \\
u(x, 0) = \sin \pi x, u_t(x, 0) = \sin \pi x
\end{cases} \tag{6}
$$

**解**：对 $t$ 作Laplace变换，得到常微分方程

$$
U_{xx} - s^2 U = -(1 + s) \sin \pi x
$$

其通解为

$$
U(x, s) = A e^{sx} + B e^{-sx} + \frac{(s + 1) \sin \pi x}{s^2 + \pi^2}
$$

由边界条件 $U(0, s) = U(1, s) = 0$ 得 $A = B = 0$，因此

$$
U(x, s) = \frac{(s + 1) \sin \pi x}{s^2 + \pi^2}
$$

作逆变换

$$
u(x, t) = \mathcal{L}^{-1}\left[\frac{(s + 1) \sin \pi x}{s^2 + \pi^2}\right] = \sin \pi x \left( \cos \pi t + \frac{1}{\pi} \sin \pi t \right)
$$

---

### 例4：半无界热传导方程

$$
\begin{cases}
u_t = u_{xx}, & x > 0, t > 0 \\
u(x, 0) = 0, & x \geq 0 \\
u(0, t) = f(t), & t > 0 \\
\lim\limits_{x \to \infty} u(x, t) = 0, & t \geq 0
\end{cases} \tag{7}
$$

**解**：设 $U(x, s) = \mathcal{L}[u(x, t)]$，则变换后方程为

$$
\frac{d^2 U}{dx^2} - s U = 0
$$

通解为

$$
U(x, s) = A e^{\sqrt{s} x} + B e^{-\sqrt{s} x}
$$

由 $x \to \infty$ 时 $U \to 0$，得 $A = 0$

由边界条件 $U(0, s) = \mathcal{L}[f(t)] = F(s)$，得 $B = F(s)$

因此

$$
U(x, s) = F(s) e^{-\sqrt{s} x}
$$

利用卷积定理作逆变换

$$
u(x, t) = \int_{0}^{t} \frac{x}{2\sqrt{\pi} \tau^{3/2}} e^{-x^2/(4\tau)} f(t - \tau) d\tau
$$

令 $\xi = x/(2\sqrt{\tau})$，则

$$
u(x, t) = \frac{2}{\sqrt{\pi}} \int_{x/(2\sqrt{t})}^{\infty} t e^{-\xi^2} f\left(t - \frac{x^2}{4\xi^2}\right) d\xi
$$

---

## 📝 考试真题与解题模板

### 真题1（2023级期末试题）

**题目**：用Laplace变换法求解波动方程半无界定解问题

$$
\begin{cases}
u_{tt} = a^2 u_{xx}, & x > 0, t > 0 \\
u(x, 0) = 0, u_t(x, 0) = 0, & x \geq 0 \\
u(0, t) = f(t), \lim\limits_{x \to +\infty} u(x, t) = 0, & t \geq 0
\end{cases}
$$

---

### 真题2（填空题常考）

求下列函数的Laplace变换：

1. $f(t) = t e^{at} \sin bt$
2. $f(t) = \int_{0}^{t} \tau e^{a\tau} \sin b\tau d\tau$

---

## 本章小结

| 性质 | 公式 | 考试频率 |
|------|------|----------|
| **微分性质** | $\mathcal{L}[f^{(n)}(t)] = s^n F(s) - s^{n-1} f(0) - \cdots - f^{(n-1)}(0)$ | ⭐⭐⭐⭐⭐ |
| **卷积定理** | $\mathcal{L}[f * g] = F(s) G(s)$ | ⭐⭐⭐⭐⭐ |
| **位移性质** | $\mathcal{L}[e^{at} f(t)] = F(s - a)$ | ⭐⭐⭐⭐ |
| **延迟性质** | $\mathcal{L}^{-1}[e^{-as} F(s)] = H(t - a) f(t - a)$ | ⭐⭐⭐⭐ |
| **乘多项式** | $\mathcal{L}[t^n f(t)] = (-1)^n F^{(n)}(s)$ | ⭐⭐⭐⭐ |

> 💡 **考试提醒**：Laplace变换每年必考！重点掌握：
> 1. 基本变换表（熟记！）
> 2. 微分性质（求解方程的核心）
> 3. 卷积定理（非齐次项处理）
> 4. 常微分方程初值问题
> 5. 半无界偏微分方程（波动/热传导）

---

> 本文件基于PPT内容整理优化，补充了详细的公式推导和解题步骤。

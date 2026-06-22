# 第4章 分离变量法

> ⭐ **考试重点**：分离变量法每年必考12分大题，特征函数法、圆域Laplace方程为高频考点

---

## 目录

- [4.1 傅里叶级数](#41-傅里叶级数)
- [4.2 齐次方程齐次边界](#42-齐次方程齐次边界)
- [4.3 非齐次方程齐次边界（特征函数法）](#43-非齐次方程齐次边界特征函数法)
- [4.4 非齐次边界处理](#44-非齐次边界处理)
- [4.5 Sturm-Liouville问题](#45-sturm-liouville问题)
- [📝 考试真题与解题模板](#-考试真题与解题模板)

---

## 4.1 傅里叶级数

### 三角函数系的正交性

三角函数系：
$$
\left\{1, \cos\frac{\pi x}{L}, \sin\frac{\pi x}{L}, \cos\frac{2\pi x}{L}, \sin\frac{2\pi x}{L}, \cdots, \cos\frac{k\pi x}{L}, \sin\frac{k\pi x}{L}, \cdots\right\}
$$

在区间 $[-L, L]$ 上是正交的。

---

### 傅里叶级数展开定理

**定理**：如果函数 $f(x)$ 以 $2L$ 为周期，在 $[-L, L]$ 上分段光滑，则对任意 $x$，成立

$$
\frac{a_0}{2} + \sum_{k=1}^{\infty} \left(a_k \cos\frac{k\pi x}{L} + b_k \sin\frac{k\pi x}{L}\right) = \frac{1}{2}\left(f(x+0) + f(x-0)\right) \tag{1}
$$

其中傅里叶系数：

$$
a_k = \frac{1}{L} \int_{-L}^{L} f(x) \cos\frac{k\pi x}{L} dx, \quad k = 0, 1, 2, \cdots \tag{2}
$$

$$
b_k = \frac{1}{L} \int_{-L}^{L} f(x) \sin\frac{k\pi x}{L} dx, \quad k = 1, 2, \cdots \tag{3}
$$

---

### 连续点的情形

如果 $f(x)$ 连续，则：
$$
f(x) = \frac{a_0}{2} + \sum_{k=1}^{\infty} \left(a_k \cos\frac{k\pi x}{L} + b_k \sin\frac{k\pi x}{L}\right) \tag{4}
$$

---

### 奇延拓（正弦级数）

如果 $f(x)$ 是奇函数（或作奇延拓），则：
$$
f(x) = \sum_{k=1}^{\infty} b_k \sin\frac{k\pi x}{L} \tag{5}
$$

其中系数：
$$
b_k = \frac{2}{L} \int_{0}^{L} f(x) \sin\frac{k\pi x}{L} dx, \quad k = 1, 2, \cdots \tag{6}
$$

> 💡 **适用：Dirichlet边界 $u(0,t)=u(L,t)=0$ 时用正弦级数

---

### 偶延拓（余弦级数）

如果 $f(x)$ 是偶函数（或作偶延拓），则：
$$
f(x) = \frac{a_0}{2} + \sum_{k=1}^{\infty} a_k \cos\frac{k\pi x}{L} \tag{7}
$$

其中系数：
$$
a_k = \frac{2}{L} \int_{0}^{L} f(x) \cos\frac{k\pi x}{L} dx, \quad k = 0, 1, 2, \cdots \tag{8}
$$

> 💡 **适用**：Neumann边界 $u_x(0,t)=u_x(L,t)=0$ 时用余弦级数

---

## 4.2 齐次方程齐次边界

### 分离变量法解题步骤（必背）

| 步骤 | 操作 |
|------|------|
| **Step 1** | 设 $u(x,t) = X(x)T(t)$，代入PDE和齐次边界条件 |
| **Step 2** | 分离变量，得到关于 $X(x)$ 的常微分方程+边界条件（特征值问题） |
| **Step 3** | 求解特征值和特征函数 |
| **Step 4** | 求解关于 $T(t)$ 的常微分方程 |
| **Step 5** | 利用叠加原理，构造级数形式解 |
| **Step 6** | 利用初始条件和正交性确定傅里叶系数 |

---

### 例1：弦振动方程（Dirichlet边界）

考虑两端固定的弦振动问题：

$$
\begin{cases}
u_{tt} = a^2 u_{xx}, & 0 < x < L, \; t > 0 \\
u(0,t) = u(L,t) = 0, & t \geq 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x \leq L
\end{cases} \tag{9}
$$

**Step 1-3. 求特征值和特征函数**

设 $u(x,t) = X(x)T(t)$，代入方程：

$$
X T'' = a^2 X'' T \implies \frac{T''}{a^2 T} = \frac{X''}{X} = -\lambda
$$

得到两个常微分方程：

$$
X'' + \lambda X = 0, \quad T'' + a^2 \lambda T = 0
$$

边界条件：$X(0) = 0, \; X(L) = 0$

**特征值**：$\lambda_n = \left(\frac{n\pi}{L}\right)^2, \quad n = 1, 2, 3, \cdots$

**特征函数**：$X_n(x) = \sin\frac{n\pi x}{L}, \quad n = 1, 2, 3, \cdots$

**Step 4. 求解时间函数**

$$
T_n(t) = A_n \cos\frac{n\pi a t}{L} + B_n \sin\frac{n\pi a t}{L}
$$

**Step 5. 叠加**

$$
u(x,t) = \sum_{n=1}^{\infty} \left(A_n \cos\frac{n\pi a t}{L} + B_n \sin\frac{n\pi a t}{L}\right) \sin\frac{n\pi x}{L}
$$

**Step 6. 确定系数**

由 $u(x,0) = \varphi(x)$：
$$
\varphi(x) = \sum_{n=1}^{\infty} A_n \sin\frac{n\pi x}{L}
$$

由 $u_t(x,0) = \psi(x)$：
$$
\psi(x) = \sum_{n=1}^{\infty} B_n \cdot \frac{n\pi a}{L} \sin\frac{n\pi x}{L}
$$

利用正交性：

$$
A_n = \frac{2}{L} \int_{0}^{L} \varphi(x) \sin\frac{n\pi x}{L} dx
$$

$$
B_n = \frac{2}{n\pi a} \int_{0}^{L} \psi(x) \sin\frac{n\pi x}{L} dx
$$

---

### 例2：弦振动方程（Neumann边界）

$$
\begin{cases}
u_{tt} = a^2 u_{xx}, & 0 < x < L, \; t > 0 \\
u_x(0,t) = u_x(L,t) = 0, & t \geq 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x \leq L
\end{cases} \tag{10}
$$

**特征值**：$\lambda_n = \left(\frac{n\pi}{L}\right)^2, \quad n = 0, 1, 2, \cdots$

**特征函数**：$X_n(x) = \cos\frac{n\pi x}{L}$

**解**：

$$
u(x,t) = \frac{A_0}{2} + \sum_{n=1}^{\infty} \left(A_n \cos\frac{n\pi a t}{L} + B_n \sin\frac{n\pi a t}{L}\right) \cos\frac{n\pi x}{L}
$$

---

### 例3：热传导方程（第三类边界）

$$
\begin{cases}
u_t = a^2 u_{xx}, & 0 < x < L, \; t > 0 \\
u(0,t) = 0, \; u_x(L,t) + h u(L,t) = 0, & t \geq 0 \\
u(x,0) = \varphi(x), & 0 \leq x \leq L
\end{cases} \tag{11}
$$

其中 $h$ 为常数。

**特征值**：$\lambda_n = \beta_n^2$，其中 $\beta_n$ 是方程 $\tan\beta L = -\frac{\beta}{h}$ 的正根

**特征函数**：$X_n(x) = \sin\beta_n x$

**解**：

$$
u(x,t) = \sum_{n=1}^{\infty} A_n e^{-a^2 \beta_n^2 t} \sin\beta_n x
$$

其中

$$
A_n = \frac{\int_{0}^{L} \varphi(x) \sin\beta_n x \, dx}{\int_{0}^{L} \sin^2\beta_n x \, dx}
$$

---

### 例4：圆域内的Dirichlet问题

$$
\begin{cases}
\Delta u = u_{xx} + u_{yy} = 0, & x^2 + y^2 < a^2 \\
u = f(x,y), & x^2 + y^2 = a^2
\end{cases} \tag{12}
$$

转换为极坐标 $(r, \theta)$：

$$
\begin{cases}
u_{rr} + \frac{1}{r} u_r + \frac{1}{r^2} u_{\theta\theta} = 0, & 0 \leq r < a \\
u(a, \theta) = f(\theta), & 0 \leq \theta \leq 2\pi
\end{cases}
$$

**分离变量**：设 $u(r, \theta) = R(r) \Theta(\theta)$

$\Theta(\theta)$ 满足：$\Theta'' + \lambda \Theta = 0$，且周期条件 $\Theta(\theta + 2\pi) = \Theta(\theta)$

$R(r)$ 满足：$r^2 R'' + r R' - \lambda R = 0$

**特征值**：$\lambda_n = n^2, \quad n = 0, 1, 2, \cdots$

**特征函数**：$\Theta_n(\theta) = a_n \cos n\theta + b_n \sin n\theta$

**径向解**：$R_0 = C_0 + D_0 \ln r$（舍去 $\ln r$ 项，因为 $r \to 0$ 时有界）
$R_n(r) = C_n r^n + D_n r^{-n}$（舍去 $r^{-n}$ 项）

**解**：

$$
u(r, \theta) = \frac{a_0}{2} + \sum_{n=1}^{\infty} \left(\frac{r}{a}\right)^n (a_n \cos n\theta + b_n \sin n\theta)
$$

其中

$$
a_n = \frac{1}{\pi} \int_{0}^{2\pi} f(\theta) \cos n\theta \, d\theta, \quad n = 0, 1, 2, \cdots
$$

$$
b_n = \frac{1}{\pi} \int_{0}^{2\pi} f(\theta) \sin n\theta \, d\theta, \quad n = 1, 2, \cdots
$$

---

## 4.3 非齐次方程齐次边界（特征函数法）

### 特征函数法（固有函数法）解题思想：

将解展开为对应齐次问题特征函数的级数

### 例：非齐次弦振动方程

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + f(x,t), & 0 < x < L, \; t > 0 \\
u(0,t) = u(L,t) = 0, & t \geq 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x \leq L
\end{cases} \tag{13}
$$

**Step 1.** 对应齐次问题的特征函数：$X_n(x) = \sin\frac{n\pi x}{L}$

**Step 2.** 设解为：
$$
u(x,t) = \sum_{n=1}^{\infty} T_n(t) \sin\frac{n\pi x}{L}
$$

**Step 3.** 将 $f(x,t)$ 也展开为特征函数级数：
$$
f(x,t) = \sum_{n=1}^{\infty} f_n(t) \sin\frac{n\pi x}{L}
$$

其中
$$
f_n(t) = \frac{2}{L} \int_{0}^{L} f(x,t) \sin\frac{n\pi x}{L} dx
$$

**Step 4.** 代入方程，得到关于 $T_n(t)$ 的常微分方程：
$$
T_n''(t) + \left(\frac{n\pi a}{L}\right)^2 T_n(t) = f_n(t)
$$

**Step 5.** 初始条件展开：
$$
\varphi(x) = \sum_{n=1}^{\infty} \varphi_n \sin\frac{n\pi x}{L}, \quad \psi(x) = \sum_{n=1}^{\infty} \psi_n \sin\frac{n\pi x}{L}
$$

得到初始条件：
$$
T_n(0) = \varphi_n, \quad T_n'(0) = \psi_n
$$

**Step 6.** 求解常微分方程初值问题，得到 $T_n(t)$

---

### 具体例题

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + \sin\frac{2\pi x}{L} \sin\frac{2\pi a t}{L}, & 0 < x < L, \; t > 0 \\
u(0,t) = u(L,t) = 0, & t \geq 0 \\
u(x,0) = 0, \; u_t(x,0) = 0, & 0 \leq x \leq L
\end{cases} \tag{14}
$$

**解**：特征函数为 $\sin\frac{n\pi x}{L}$

设
$$
u(x,t) = \sum_{n=1}^{\infty} T_n(t) \sin\frac{n\pi x}{L}
$$

非齐次项已是特征函数，仅 $n=2$ 项非零：

$$
T_2''(t) + \left(\frac{2\pi a}{L}\right)^2 T_2(t) = \sin\frac{2\pi a t}{L}
$$

初始条件：$T_2(0) = 0, \; T_2'(0) = 0$

用常数变易法或拉普拉斯变换求解：

$$
T_2(t) = \frac{L}{4\pi a} \left(\frac{L}{2\pi a} \sin\frac{2\pi a t}{L} - t \cos\frac{2\pi a t}{L}\right)
$$

最终解：

$$
u(x,t) = \frac{L^2}{8\pi^2 a^2} \left(\sin\frac{2\pi a t}{L} - \frac{2\pi a t}{L} \cos\frac{2\pi a t}{L}\right) \sin\frac{2\pi x}{L}
$$

---

## 4.4 非齐次边界处理

### 基本思想：通过变量替换将非齐次边界齐次化

考虑问题：

$$
\begin{cases}
u_{tt} = a^2 u_{xx} + f(x,t), & 0 < x < L, \; t > 0 \\
u(0,t) = p(t), \; u(L,t) = q(t), & t \geq 0 \\
u(x,0) = \varphi(x), \; u_t(x,0) = \psi(x), & 0 \leq x \leq L
\end{cases} \tag{15}
$$

**Step 1.** 构造辅助函数 $v(x,t)$ 满足非齐次边界条件，通常取线性函数：

$$
v(x,t) = p(t) + \frac{x}{L} \left(q(t) - p(t)\right)
$$

**Step 2.** 令 $u(x,t) = v(x,t) + w(x,t)$，则 $w(x,t)$ 满足齐次边界：

$$
w(0,t) = w(L,t) = 0
$$

**Step 3.** 推导 $w$ 满足的方程和初始条件，用特征函数法求解

---

## 4.5 Sturm-Liouville问题

### 一般形式

$$
\begin{cases}
\frac{d}{dx} \left[p(x) X'(x)\right] + \left[q(x) + \lambda r(x)\right] X(x) = 0, \quad a < x < b \\
\alpha_1 X(a) + \alpha_2 X'(a) = 0 \\
\beta_1 X(b) + \beta_2 X'(b) = 0
\end{cases}
$$

其中 $p(x) > 0, \; r(x) > 0$

### 重要性质（Sturm-Liouville定理）

1. **特征值存在且都是实数**，可排列为：$\lambda_1 < \lambda_2 < \cdots < \lambda_n < \cdots \to +\infty$
2. **特征函数系 $\{X_n(x)\}$ 关于权函数 $r(x)$ 正交**：
   $$
   \int_{a}^{b} r(x) X_m(x) X_n(x) dx = 0, \quad m \neq n
   $$
3. **特征函数系在加权空间中是完备的**，任意平方可积函数可展开为特征函数的傅里叶级数

> 💡 **分离变量法的理论基础**：保证解的存在性、正交性、完备性

---

## 📝 考试真题与解题模板

### 分离变量法解题步骤（必背）

| 问题类型 | 边界条件 | 特征函数 | 解法要点 |
|----------|----------|----------|----------|
| 波动方程 | $u(0)=u(L)=0$ | $\sin\frac{n\pi x}{L}$ | 时间函数含 $\cos, \sin$ |
| 热传导方程 | $u(0)=u(L)=0$ | $\sin\frac{n\pi x}{L}$ | 时间函数含 $e^{-a^2\lambda t}$ |
| Laplace方程 | 圆域边界 | $r^n \cos n\theta, r^n \sin n\theta$ | 舍去无界项 |

---

### 真题1（2023级期末试题）⭐⭐⭐⭐⭐

**题目**：用特征函数法求解热传导方程混合问题

$$
\begin{cases}
u_t = u_{xx} + \sin\frac{x}{2}, & 0 < x < \pi, \; t > 0 \\
u(0,t) = 0, \; u_x(\pi,t) = 0, & t \geq 0 \\
u(x,0) = \sin\frac{x}{2}, & 0 \leq x \leq \pi
\end{cases}
$$

**解答：**

**Step 1. 求对应齐次问题的特征值和特征函数**

齐次问题：
$$
\begin{cases}
X'' + \lambda X = 0 \\
X(0) = 0, \; X'(\pi) = 0
\end{cases}
$$

通解：$X(x) = A\cos\sqrt{\lambda}x + B\sin\sqrt{\lambda}x$

由 $X(0) = 0 \implies A = 0$

由 $X'(\pi) = B\sqrt{\lambda}\cos\sqrt{\lambda}\pi = 0 \implies \cos\sqrt{\lambda}\pi = 0$

$$
\sqrt{\lambda}\pi = \frac{(2n-1)\pi}{2} \implies \lambda_n = \left(\frac{2n-1}{2}\right)^2, \quad n = 1, 2, \cdots
$$

特征函数：
$$
X_n(x) = \sin\frac{(2n-1)x}{2}
$$

**Step 2. 设解为特征函数级数**

$$
u(x,t) = \sum_{n=1}^{\infty} T_n(t) \sin\frac{(2n-1)x}{2}
$$

**Step 3. 将非齐次项展开**

注意到 $\sin\frac{x}{2}$ 本身就是 $n=1$ 的特征函数

**Step 4. 代入方程**

$$
T_n'(t) + \lambda_n T_n(t) = 
\begin{cases}
1, & n = 1 \\
0, & n \geq 2
\end{cases}
$$

**Step 5. 初始条件**

$$
u(x,0) = \sin\frac{x}{2} = \sum_{n=1}^{\infty} T_n(0) \sin\frac{(2n-1)x}{2}
$$

所以 $T_1(0) = 1, \; T_n(0) = 0 \; (n \geq 2)$

**Step 6. 求解常微分方程**

对 $n = 1$：
$$
T_1'(t) + \frac{1}{4} T_1(t) = 1, \quad T_1(0) = 1
$$

解得：
$$
T_1(t) = 4 - 3e^{-t/4}
$$

对 $n \geq 2$：
$$
T_n'(t) + \lambda_n T_n(t) = 0, \quad T_n(0) = 0 \implies T_n(t) = 0
$$

**Step 7. 最终解**

$$
u(x,t) = (4 - 3e^{-t/4}) \sin\frac{x}{2}
$$

---

### 真题2（扇形区域Laplace方程）⭐⭐⭐⭐

**题目**：求解扇形区域上拉普拉斯方程边值问题

$$
\begin{cases}
\Delta u = u_{rr} + \frac{1}{r} u_r + \frac{1}{r^2} u_{\theta\theta} = 0, & r < 1, \; 0 < \theta < \frac{\pi}{4} \\
u(r,0) = 0, \quad u_\theta\left(r, \frac{\pi}{4}\right) = 0, & 0 \leq r \leq 1 \\
u(1,\theta) = \sin 2\theta + \sin 6\theta, & 0 \leq \theta \leq \frac{\pi}{4}
\end{cases}
$$

**解答：**

**Step 1. 分离变量**

设 $u(r, \theta) = R(r) \Theta(\theta)$，代入方程得：

$$
\frac{r^2 R'' + r R'}{R} = -\frac{\Theta''}{\Theta} = \lambda
$$

得到两个常微分方程：

$$
r^2 R'' + r R' - \lambda R = 0
$$
$$
\Theta'' + \lambda \Theta = 0
$$

**Step 2. 求解特征值问题**

边界条件：$\Theta(0) = 0, \; \Theta'\left(\frac{\pi}{4}\right) = 0$

通解：$\Theta(\theta) = A\cos\sqrt{\lambda}\theta + B\sin\sqrt{\lambda}\theta$

由 $\Theta(0) = 0 \implies A = 0$

由 $\Theta'\left(\frac{\pi}{4}\right) = B\sqrt{\lambda}\cos\left(\sqrt{\lambda} \cdot \frac{\pi}{4}\right) = 0$

$$
\sqrt{\lambda} \cdot \frac{\pi}{4} = \frac{(2n-1)\pi}{2} \implies \lambda_n = 2(2n-1))^2, \quad n = 1, 2, \cdots
$$

特征函数：
$$
\Theta_n(\theta) = \sin 2(2n-1)\theta
$$

**Step 3. 求解径向方程**

欧拉方程解为：
$$
R_n(r) = C_n r^{2(2n-1)} + D_n r^{-2(2n-1)}
$$

舍去 $r^{-2(2n-1)}$ 项（$r \to 0$ 时无界），故
$$
R_n(r) = C_n r^{2(2n-1)}
$$

**Step 4. 叠加**

$$
u(r, \theta) = \sum_{n=1}^{\infty} C_n r^{2(2n-1)} \sin 2(2n-1)\theta
$$

**Step 5. 利用边界条件确定系数**

$$
u(1, \theta) = \sin 2\theta + \sin 6\theta = \sum_{n=1}^{\infty} C_n \sin 2(2n-1)\theta
$$

对比得：
- $n=1$ 时：$C_1 = 1$（对应 $\sin 2\theta$）
- $n=2$ 时：$C_2 = 1$（对应 $\sin 6\theta$）
- 其他 $C_n = 0$

**Step 6. 最终解**

$$
u(r, \theta) = r^2 \sin 2\theta + r^6 \sin 6\theta
$$

---

## 本章小结

| 方法 | 适用问题 | 核心思想 | 考试频率 |
|------|----------|----------|----------|
| **分离变量法** | 齐次方程齐次边界 | 设u=X(x)T(t)，解特征值问题 | ⭐⭐⭐⭐⭐ |
| **特征函数法** | 非齐次方程齐次边界 | 解展开为齐次问题特征函数级数 | ⭐⭐⭐⭐⭐ |
| **边界齐次化** | 非齐次边界 | 构造辅助函数消去非齐次性 | ⭐⭐⭐ |
| **极坐标分离变量 | 圆域、扇形Laplace方程 | 欧拉方程+周期/边界条件 | ⭐⭐⭐⭐ |

> 💡 **考试提醒**：分离变量法每年必考12分大题！重点掌握：(1) 热传导方程特征函数法（2023级已考）；(2) 圆域/扇形Laplace方程（高频）；(3) 波动方程分离变量。

---

> 本文件基于PPT内容整理优化，补充了详细的解题步骤和历年考试真题详解。

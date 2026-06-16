# 题型 C 详细解法解释：圆域 / 极坐标 Laplace / Poisson 方程

## 一、为什么需要极坐标？

当区域是**圆形**时（如圆盘 $r < R$），用直角坐标 $(x, y)$ 很难处理边界条件，因为圆的边界 $x^2 + y^2 = R^2$ 不是坐标线。

所以需要换成**极坐标**：
$$x = r\cos\theta, \quad y = r\sin\theta$$

在极坐标下，圆的边界就是 $r = R$（常数），处理起来方便多了。

---

## 二、极坐标下的拉普拉斯算子

在直角坐标中：
$$\Delta u = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}$$

在极坐标中（需要推导）：
$$\boxed{\Delta u = \frac{\partial^2 u}{\partial r^2} + \frac{1}{r}\frac{\partial u}{\partial r} + \frac{1}{r^2}\frac{\partial^2 u}{\partial \theta^2}}$$

> **推导思路**：用链式法则，把 $u_{xx} + u_{yy}$ 用 $r, \theta$ 表示。这里不展开推导，直接记住结论。

---

## 三、极坐标下的分离变量

### 3.1 设分离变量形式

对于 Laplace 方程 $\Delta u = 0$，设：
$$u(r, \theta) = R(r) \cdot \Theta(\theta)$$

### 3.2 代入拉普拉斯算子

$$\frac{\partial^2 u}{\partial r^2} = R''(r)\Theta(\theta)$$

$$\frac{1}{r}\frac{\partial u}{\partial r} = \frac{1}{r}R'(r)\Theta(\theta)$$

$$\frac{1}{r^2}\frac{\partial^2 u}{\partial \theta^2} = \frac{1}{r^2}R(r)\Theta''(\theta)$$

代入 $\Delta u = 0$：
$$R''\Theta + \frac{1}{r}R'\Theta + \frac{1}{r^2}R\Theta'' = 0$$

### 3.3 分离变量

两边乘以 $\frac{r^2}{R\Theta}$：
$$\frac{r^2 R'' + rR'}{R} = -\frac{\Theta''}{\Theta} = \lambda$$

得到两个常微分方程：
1. **角向方程**：$\Theta'' + \lambda \Theta = 0$
2. **径向方程**：$r^2 R'' + rR' - \lambda R = 0$（欧拉方程）

---

## 四、角向方程的求解（关键！）

### 4.1 周期性边界条件

角向方程 $\Theta'' + \lambda \Theta = 0$ 需要满足**周期性条件**：
$$\Theta(\theta) = \Theta(\theta + 2\pi)$$

因为 $\theta$ 和 $\theta + 2\pi$ 代表同一个方向。

### 4.2 求解特征值

设 $\lambda > 0$，令 $\lambda = n^2$：
$$\Theta'' + n^2 \Theta = 0$$

通解：
$$\Theta(\theta) = A\cos(n\theta) + B\sin(n\theta)$$

**周期性要求**：$n$ 必须是**正整数**（$n = 1, 2, 3, \ldots$），否则 $\cos(n\theta)$ 不是周期的。

当 $\lambda = 0$ 时：$\Theta = A + B\theta$，周期性要求 $B = 0$，所以 $\Theta = 1$（常数）。

### 4.3 特征值和特征函数

| 特征值 $\lambda$ | 特征函数 $\Theta(\theta)$ |
|-----------------|-------------------------|
| $\lambda_0 = 0$ | $\Theta_0 = 1$ |
| $\lambda_1 = 1^2$ | $\cos\theta, \sin\theta$ |
| $\lambda_2 = 2^2$ | $\cos 2\theta, \sin 2\theta$ |
| $\lambda_n = n^2$ | $\cos(n\theta), \sin(n\theta)$ |

---

## 五、径向方程的求解

### 5.1 欧拉方程

径向方程 $r^2 R'' + rR' - n^2 R = 0$ 是**欧拉方程**，设 $R = r^m$：
$$m(m-1)r^m + mr^m - n^2 r^m = 0$$
$$m^2 - n^2 = 0 \implies m = \pm n$$

### 5.2 通解形式

- 当 $n = 0$：$R = C_0 + D_0 \ln r$
- 当 $n \geq 1$：$R = C_n r^n + D_n r^{-n}$

### 5.3 正则性条件（重要！）

圆心 $r = 0$ 必须在区域内，所以解在 $r = 0$ 处必须**有界**：
- 舍去 $\ln r$（$n=0$ 时）
- 舍去 $r^{-n}$（$n \geq 1$ 时）

因此：
- $n = 0$：$R_0 = 1$（常数）
- $n \geq 1$：$R_n = r^n$

---

## 六、调和函数的一般形式（核心结论！）

把角向和径向的解组合起来：

$$\boxed{u(r, \theta) = \frac{a_0}{2} + \sum_{n=1}^{\infty} r^n \left[a_n \cos(n\theta) + b_n \sin(n\theta)\right]}$$

这就是**圆盘内调和函数**（$\Delta u = 0$）的一般形式。

> **注意**：
> - $\frac{a_0}{2}$ 中的 $\frac{1}{2}$ 是为了后面傅里叶系数公式简洁
> - $r^n$ 保证在 $r = 0$ 处有界

---

## 七、应用边界条件

### 7.1 边界条件

在圆的边界 $r = R$ 处，给定：
$$u(R, \theta) = f(\theta)$$

### 7.2 代入一般形式

$$f(\theta) = \frac{a_0}{2} + \sum_{n=1}^{\infty} R^n \left[a_n \cos(n\theta) + b_n \sin(n\theta)\right]$$

### 7.3 傅里叶系数公式

这就是 $f(\theta)$ 的**傅里叶级数展开**，系数为：

$$a_0 = \frac{1}{\pi} \int_0^{2\pi} f(\theta) \, d\theta$$

$$a_n = \frac{1}{\pi R^n} \int_0^{2\pi} f(\theta) \cos(n\theta) \, d\theta$$

$$b_n = \frac{1}{\pi R^n} \int_0^{2\pi} f(\theta) \sin(n\theta) \, d\theta$$

---

## 八、Poisson 方程的处理

对于 $\Delta u = f$（非齐次方程），需要用**特解 + 调和函数**的方法。

### 8.1 为什么要这样做？

- $\Delta u = f$ 不能直接用分离变量（因为有源项 $f$）
- 先找一个**特解 $u_p$** 满足 $\Delta u_p = f$（不管边界条件）
- 然后 $v = u - u_p$ 满足 $\Delta v = 0$（调和函数，用分离变量）

### 8.2 特解的构造（公式法）

利用这个**神奇公式**：
$$\boxed{\Delta(r^k \cos m\theta) = (k^2 - m^2) r^{k-2} \cos m\theta}$$

以及：
$$\Delta(r^k \sin m\theta) = (k^2 - m^2) r^{k-2} \sin m\theta$$

**这个公式的推导**（帮助理解）：

在极坐标下计算：
$$\Delta(r^k \cos m\theta) = \frac{\partial^2}{\partial r^2}(r^k \cos m\theta) + \frac{1}{r}\frac{\partial}{\partial r}(r^k \cos m\theta) + \frac{1}{r^2}\frac{\partial^2}{\partial \theta^2}(r^k \cos m\theta)$$

$$= k(k-1)r^{k-2}\cos m\theta + kr^{k-2}\cos m\theta + \frac{1}{r^2}(-m^2 r^k \cos m\theta)$$

$$= [k(k-1) + k - m^2] r^{k-2} \cos m\theta = (k^2 - m^2) r^{k-2} \cos m\theta$$

### 8.3 如何用公式构造特解？

**例子**：源项 $f = r^2$

观察：$r^2 = r^2 \cos(0\theta)$，即 $k = 2, m = 0$。

$$\Delta(r^2 \cos 0\theta) = (2^2 - 0^2) r^{2-2} \cos 0\theta = 4$$

所以 $\Delta\left(\frac{r^2}{4}\right) = 1$。

因此，$\Delta u = r^2$ 的一个特解是 $u_p = \frac{r^4}{16}$。

验证：$\Delta(r^4) = (4^2 - 0) r^2 = 16 r^2$，所以 $\Delta\left(\frac{r^4}{16}\right) = r^2$ ✓

---

## 九、完整例题：一步步详解

### 例题

求解 Poisson 方程：
$$\begin{cases}
\Delta u = r^2, & r < 1 \\
u(1, \theta) = \cos(2\theta)
\end{cases}$$

### Step 1：找特解 $u_p$（满足方程，不管边界）

源项是 $r^2 = r^2 \cos(0\theta)$，对应 $k=2, m=0$。

用公式：$\Delta(r^4) = 16 r^2$，所以：
$$u_p = \frac{r^4}{16}$$

验证：$\Delta u_p = \Delta\left(\frac{r^4}{16}\right) = \frac{16r^2}{16} = r^2$ ✓

### Step 2：找调和函数 $v = u - u_p$

令 $v = u - u_p$，则：
$$\Delta v = \Delta u - \Delta u_p = r^2 - r^2 = 0$$

所以 $v$ 是调和函数！

### Step 3：写出 $v$ 的一般形式

$$v(r, \theta) = \frac{a_0}{2} + \sum_{n=1}^{\infty} r^n [a_n \cos(n\theta) + b_n \sin(n\theta)]$$

### Step 4：应用边界条件

在 $r = 1$ 处：
$$u(1, \theta) = u_p(1) + v(1, \theta)$$

$$\cos(2\theta) = \frac{1}{16} + v(1, \theta)$$

所以：
$$v(1, \theta) = \cos(2\theta) - \frac{1}{16}$$

### Step 5：确定 $v$ 的系数

在 $r = 1$ 处展开：
$$v(1, \theta) = \frac{a_0}{2} + \sum_{n=1}^{\infty} [a_n \cos(n\theta) + b_n \sin(n\theta)]$$

比较系数：
- 常数项：$\frac{a_0}{2} = -\frac{1}{16} \implies a_0 = -\frac{1}{8}$
- $n=2, \cos$ 项：$a_2 = 1$
- 其他所有系数：$0$

### Step 6：写出 $v$

$$v(r, \theta) = -\frac{1}{16} + r^2 \cos(2\theta)$$

### Step 7：合成最终解

$$u = u_p + v = \frac{r^4}{16} - \frac{1}{16} + r^2 \cos(2\theta)$$

$$\boxed{u(r, \theta) = \frac{r^4 - 1}{16} + r^2 \cos(2\theta)}$$

### Step 8：验证

1. **边界条件**：$r=1$ 时，$u = \frac{1-1}{16} + 1 \cdot \cos(2\theta) = \cos(2\theta)$ ✓

2. **方程**：
   - $u_p = \frac{r^4}{16}$，$\Delta u_p = r^2$
   - $v = -\frac{1}{16} + r^2 \cos(2\theta)$
   - $\Delta v = \Delta(r^2 \cos 2\theta) = (4-4) r^0 \cos 2\theta = 0$ ✓
   - 所以 $\Delta u = r^2$ ✓

---

## 十、总结：题型 C 的解题步骤

### 对于纯 Laplace 方程 $\Delta u = 0$

1. **写出调和函数的一般形式**：
   $$u(r, \theta) = \frac{a_0}{2} + \sum_{n=1}^{\infty} r^n [a_n \cos(n\theta) + b_n \sin(n\theta)]$$

2. **代入边界条件 $u(R, \theta) = f(\theta)$**

3. **计算傅里叶系数**：
   $$a_n = \frac{1}{\pi R^n} \int_0^{2\pi} f(\theta) \cos(n\theta) \, d\theta$$
   $$b_n = \frac{1}{\pi R^n} \int_0^{2\pi} f(\theta) \sin(n\theta) \, d\theta$$

4. **写出最终解**

### 对于 Poisson 方程 $\Delta u = f$

1. **找特解 $u_p$**：用公式 $\Delta(r^k \cos m\theta) = (k^2-m^2)r^{k-2}\cos m\theta$ 构造

2. **令 $v = u - u_p$**：$v$ 满足 $\Delta v = 0$

3. **写出 $v$ 的一般形式**（调和函数）

4. **应用边界条件**：在 $r = R$ 处匹配系数

5. **合成最终解**：$u = u_p + v$

---

## 十一、和直角坐标分离变量的对比

| 步骤 | 直角坐标 | 极坐标 |
|------|---------|--------|
| 分离变量 | $u = X(x)T(t)$ | $u = R(r)\Theta(\theta)$ |
| 特征值问题 | $X'' + \lambda X = 0$ | $\Theta'' + n^2 \Theta = 0$ |
| 边界条件 | 两端固定 | 周期性 |
| 特征函数 | $\sin\frac{n\pi x}{L}$ | $\cos(n\theta), \sin(n\theta)$ |
| 径向/时间方程 | 常系数 ODE | 欧拉方程 |
| 径向/时间解 | 指数/三角函数 | $r^n$（幂函数） |
| 展开什么 | 初始条件 | 边界条件 |

**核心区别**：
- 直角坐标：在**初始时刻**展开
- 极坐标：在**边界**展开

---

## 十二、常见问题

### Q1：为什么 $r^{-n}$ 被舍去？

因为圆心 $r = 0$ 在区域内，$r^{-n}$ 在 $r = 0$ 处发散（无穷大），不满足物理上的有界性要求。

### Q2：为什么 $n=0$ 时是常数而不是 $\ln r$？

$\ln r$ 在 $r = 0$ 处也发散，所以舍去。$n=0$ 对应的径向方程是 $r^2 R'' + rR' = 0$，解为 $R = C_0 + D_0 \ln r$，舍去 $\ln r$ 后只剩常数 $C_0$。

### Q3：特解不唯一怎么办？

特解的选择不唯一，但最终结果 $u = u_p + v$ 是唯一的。选择最简单的特解即可。

### Q4：如何快速判断特解的形式？

观察源项 $f$ 的结构：
- 常数 $C$：$u_p = \frac{C}{4}r^2$
- $r^k \cos(m\theta)$：$u_p = \frac{r^{k+2}}{(k+2)^2-m^2} \cos(m\theta)$（当 $(k+2)^2 \neq m^2$）

---

## 十三、公式速查表

| 源项 $f$ | 特解 $u_p$ | 验证 |
|---------|-----------|------|
| 常数 $C$ | $\frac{C}{4}r^2$ | $\Delta(r^2) = 4$ |
| $r^2$ | $\frac{r^4}{16}$ | $\Delta(r^4) = 16r^2$ |
| $r^2 \cos(2\theta)$ | $\frac{r^4}{12}\cos(2\theta)$ | $\Delta(r^4\cos 2\theta) = 12r^2\cos 2\theta$ |
| $r\cos\theta$ | $\frac{r^3}{8}\cos\theta$ | $\Delta(r^3\cos\theta) = 8r\cos\theta$ |
| $\cos(2\theta)$ | $\frac{r^2}{4}\cos(2\theta)$ | $\Delta(r^2\cos 2\theta) = 0$（已是调和函数） |

---

## 十四、扇形区域问题（重要变体！）

### 14.1 与完整圆盘的区别

| 区域 | 角向范围 | 边界条件类型 | 特征值 |
|------|---------|-------------|--------|
| **完整圆盘** | $0 \leq \theta < 2\pi$ | 周期性：$\Theta(\theta) = \Theta(\theta + 2\pi)$ | $\lambda_n = n^2$，$n = 0, 1, 2, \ldots$ |
| **扇形区域** | $0 < \theta < \alpha$ | 两端固定：$\Theta(0) = 0$，$\Theta'(\alpha) = 0$（或反之） | $\lambda_n = \left(\frac{(2n-1)\pi}{2\alpha}\right)^2$ |

### 14.2 扇形区域的角向特征值问题

**典型边界条件**：
$$\Theta(0) = 0, \quad \Theta'(\alpha) = 0$$

**求解**：

设 $\lambda > 0$，令 $\lambda = \mu^2$：
$$\Theta = A\cos(\mu\theta) + B\sin(\mu\theta)$$

$\Theta(0) = 0 \implies A = 0$，所以 $\Theta = B\sin(\mu\theta)$

$\Theta'(\alpha) = B\mu\cos(\mu\alpha) = 0 \implies \cos(\mu\alpha) = 0$

$$\mu\alpha = \frac{(2n-1)\pi}{2} \implies \mu_n = \frac{(2n-1)\pi}{2\alpha}$$

**特征值和特征函数**：
$$\boxed{\lambda_n = \left(\frac{(2n-1)\pi}{2\alpha}\right)^2, \quad \Theta_n(\theta) = \sin\left(\frac{(2n-1)\pi\theta}{2\alpha}\right)}$$

### 14.3 常见扇形角度的特征值

| 扇形角度 $\alpha$ | 特征值 $\lambda_n$ | 特征函数 $\Theta_n$ |
|------------------|-------------------|---------------------|
| $\frac{\pi}{4}$ | $4(2n-1)^2$ | $\sin[2(2n-1)\theta]$ |
| $\frac{\pi}{2}$ | $(2n-1)^2$ | $\sin[(2n-1)\theta]$ |
| $\pi$ | $\frac{(2n-1)^2}{4}$ | $\sin\frac{(2n-1)\theta}{2}$ |

### 14.4 完整例题：2023级期末第3题

**题目**：
$$\begin{cases}
\Delta u = 0, & r < 1, \ 0 < \theta < \dfrac{\pi}{4} \\
u(r,0) = 0, \ u_\theta\left(r,\dfrac{\pi}{4}\right) = 0 \\
u(1,\theta) = \sin 2\theta + \sin 6\theta
\end{cases}$$

**解**：

**Step 1：角向特征值问题**

$\Theta(0) = 0$，$\Theta'(\frac{\pi}{4}) = 0$，$\alpha = \frac{\pi}{4}$

$$\mu_n = \frac{(2n-1)\pi}{2 \cdot \frac{\pi}{4}} = 2(2n-1)$$

$$\lambda_n = 4(2n-1)^2, \quad \Theta_n = \sin[2(2n-1)\theta]$$

验证：$\sin[2(2\cdot1-1)\theta] = \sin 2\theta$，$\sin[2(2\cdot2-1)\theta] = \sin 6\theta$ ✓

**Step 2：径向方程**

$$r^2R'' + rR' - 4(2n-1)^2 R = 0$$

欧拉方程：$m^2 - 4(2n-1)^2 = 0 \implies m = \pm 2(2n-1)$

舍去负幂项：$R_n = r^{2(2n-1)}$

**Step 3：一般解**

$$u = \sum_{n=1}^{\infty} A_n r^{2(2n-1)} \sin[2(2n-1)\theta]$$

**Step 4：确定系数**

$u(1,\theta) = \sin 2\theta + \sin 6\theta$

- $n=1$：$\sin 2\theta$ 的系数 $A_1 = 1$
- $n=2$：$\sin 6\theta$ 的系数 $A_2 = 1$

**最终解**：
$$\boxed{u(r,\theta) = r^2\sin 2\theta + r^6\sin 6\theta}$$

### 14.5 其他边界条件组合

| 边界条件 | 特征值 $\mu_n$ | 特征函数 |
|---------|---------------|---------|
| $\Theta(0) = 0$，$\Theta'(\alpha) = 0$ | $\frac{(2n-1)\pi}{2\alpha}$ | $\sin(\mu_n\theta)$ |
| $\Theta'(0) = 0$，$\Theta(\alpha) = 0$ | $\frac{(2n-1)\pi}{2\alpha}$ | $\cos(\mu_n\theta)$ |
| $\Theta(0) = 0$，$\Theta(\alpha) = 0$ | $\frac{n\pi}{\alpha}$ | $\sin(\mu_n\theta)$ |
| $\Theta'(0) = 0$，$\Theta'(\alpha) = 0$ | $\frac{n\pi}{\alpha}$（含 $n=0$） | $\cos(\mu_n\theta)$ |

---

## 十五、圆域问题解法总结流程图

```
开始
  ↓
判断区域类型
  ↓
┌─────────────────────────────────────┐
│ 完整圆盘（0 ≤ θ < 2π）              │
│ 角向：周期性条件 → λ = n²           │
│ 径向：R = rⁿ（舍去 r⁻ⁿ）           │
└─────────────────────────────────────┘
  ↓
┌─────────────────────────────────────┐
│ 扇形区域（0 < θ < α）               │
│ 角向：两端条件 → μ = (2n-1)π/(2α)  │
│ 径向：R = rᵐ（舍去 r⁻ᵐ）           │
└─────────────────────────────────────┘
  ↓
写出一般解
  ↓
应用边界条件 r = R 处的条件
  ↓
确定傅里叶系数
  ↓
写出最终解
```

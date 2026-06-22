# 第2章 二阶线性PDE分类与标准形

> ⭐ **考试重点**：化标准型为必考解答题（12分），三类方程的典型代表为高频填空题

---

## 目录

- [2.1 两个自变量的二阶线性偏微分方程的分类](#21-两个自变量的二阶线性偏微分方程的分类)
- [2.2 多个自变量的二阶线性偏微分方程的分类](#22-多个自变量的二阶线性偏微分方程的分类)
- [📝 考试真题与解题模板](#-考试真题与解题模板)

---

## 2.1 两个自变量的二阶线性偏微分方程的分类

### 二阶线性偏微分方程的一般形式

两个自变量的二阶线性偏微分方程的一般形式为：

$$
L[u] = a_{11}u_{xx} + 2a_{12}u_{xy} + a_{22}u_{yy} + a_1u_x + b_1u_y + cu = f \tag{1}
$$

其中 $a_{ij}, a_1, b_1, c$ 是 $(x,y)$ 的已知函数，$a_{ij}$ 不全为零 $(i,j=1,2)$。

---

### 自变量变换

作可逆的自变量变换：

$$
\xi = \xi(x,y), \quad \eta = \eta(x,y), \quad |J| = \left| \frac{\partial(\xi,\eta)}{\partial(x,y)} \right| \neq 0
$$

在此变换下，$u=u(x,y)$ 视为 $u=u(\xi,\eta)$ 的函数，方程 (1) 变为：

$$
A_{11}u_{\xi\xi} + 2A_{12}u_{\xi\eta} + A_{22}u_{\eta\eta} + A_1u_\xi + B_1u_\eta + cu = f \tag{2}
$$

其中系数的变换关系为：

$$
A_{11} = a_{11}\xi_x^2 + 2a_{12}\xi_x\xi_y + a_{22}\xi_y^2 \tag{3}
$$

$$
A_{12} = a_{11}\xi_x\eta_x + a_{12}(\xi_x\eta_y + \xi_y\eta_x) + a_{22}\xi_y\eta_y \tag{4}
$$

$$
A_{22} = a_{11}\eta_x^2 + 2a_{12}\eta_x\eta_y + a_{22}\eta_y^2 \tag{5}
$$

$$
A_1 = L[\xi] - c\xi, \quad B_1 = L[\eta] - c\eta
$$

---

### 方程的分类 ✨ 考试重点（填空题）

引入判别式：

$$
\Delta = a_{12}^2(x,y) - a_{11}(x,y)a_{22}(x,y)
$$

根据判别式的符号，在点 $(x_0,y_0)$ 处方程可分为三类：

| 类型 | 判别式条件 | 典型方程 |
|------|------------|----------|
| **双曲型** | $\Delta > 0$ | **波动方程** $u_{tt} = a^2u_{xx}$ |
| **抛物型** | $\Delta = 0$ | **热传导方程** $u_t = a^2u_{xx}$ |
| **椭圆型** | $\Delta < 0$ | **位势方程(Laplace)** $u_{xx} + u_{yy} = 0$ |

> 💡 **记忆技巧**：**波**浪有两个方向（双曲），**热**量单向传播（抛物），**位势**是封闭曲线（椭圆）

---

### 特征方程

由 (3) 和 (5) 式可见，若要使 $A_{11}$ 或 $A_{22}$ 等于零，则 $\xi$ 或 $\eta$ 必须是下列一阶偏微分方程的解：

$$
a_{11}\varphi_x^2 + 2a_{12}\varphi_x\varphi_y + a_{22}\varphi_y^2 = 0 \tag{6}
$$

这个方程等价于常微分方程：

$$
a_{11}(dy)^2 - 2a_{12}dxdy + a_{22}(dx)^2 = 0 \tag{7}
$$

称为方程 (1) 的 **特征方程**，由它所确定的方向 $\frac{dy}{dx}$ 称为 **特征方向**，其积分曲线称为方程 (1) 的 **特征曲线**。

求解特征方程可得两个微分方程：

$$
\frac{dy}{dx} = \frac{a_{12} + \sqrt{a_{12}^2 - a_{11}a_{22}}}{a_{11}} \tag{8}
$$

$$
\frac{dy}{dx} = \frac{a_{12} - \sqrt{a_{12}^2 - a_{11}a_{22}}}{a_{11}} \tag{9}
$$

---

### 三类方程的标准形式

#### 1. 双曲型方程（$\Delta > 0$）

当 $\Delta > 0$ 时，方程 (1) 是双曲型的。特征方程 (7) 有两个不同的特征方向：

$$
\frac{dy}{dx} = \frac{a_{12} + \sqrt{a_{12}^2 - a_{11}a_{22}}}{a_{11}} \tag{10}
$$

$$
\frac{dy}{dx} = \frac{a_{12} - \sqrt{a_{12}^2 - a_{11}a_{22}}}{a_{11}} \tag{11}
$$

得到两族不同的特征曲线：$\varphi(x,y) = c_1, \ \psi(x,y) = c_2$

令 $\xi = \varphi(x,y), \ \eta = \psi(x,y)$，则方程 (1) 可化为**第一标准形**：

$$
u_{\xi\eta} + A_2u_\xi + B_2u_\eta + C_2u + F_1 = 0 \tag{12}
$$

再作变换 $s = \xi + \eta, \ t = \xi - \eta$，则 (12) 式可化为**第二标准形**：

$$
u_{ss} - u_{tt} + A_3u_s + B_3u_t + C_3u + F_2 = 0 \tag{13}
$$

---

#### 2. 抛物型方程（$\Delta = 0$）

当 $\Delta = 0$ 时，方程 (1) 是抛物型的。此时只有一个特征方向：

$$
\frac{dy}{dx} = \frac{a_{12}}{a_{11}} \tag{14}
$$

特征曲线为 $\varphi(x,y) = c_1$。

令 $\xi = \varphi(x,y)$，再任意选取 $\eta = \eta(x,y)$，使得 $\left| \frac{\partial(\xi,\eta)}{\partial(x,y)} \right| \neq 0$。通常可取 $\eta = x$。

此时方程可化为：

$$
u_{\eta\eta} + A_4u_\xi + B_4u_\eta + C_4u + F_3 = 0 \tag{15}
$$

这就是抛物型方程的标准形。

---

#### 3. 椭圆型方程（$\Delta < 0$）

当 $\Delta < 0$ 时，方程 (1) 是椭圆型的。特征方程 (7) 不存在实的特征曲线，方程 (8)-(9) 的积分曲线是一对共轭的复特征曲线：

$$
\varphi(x,y) = \alpha + i\beta = c_1, \quad \psi(x,y) = \alpha - i\beta = c_2
$$

其中 $\alpha, \beta$ 是 $(x,y)$ 的实函数。

令 $\xi = \alpha, \ \eta = \beta$，则方程 (1) 可以化为：

$$
u_{\alpha\alpha} + u_{\beta\beta} + A_5u_\alpha + B_5u_\beta + C_5u + F_4 = 0 \tag{16}
$$

这就是椭圆型方程的标准形。

---

## 2.2 多个自变量的二阶线性偏微分方程的分类

> 注：本PPT该节仅在目录中列出，具体内容未在7页PPT中展开。考试主要考察两个自变量的情况。

---

## 📝 考试真题与解题模板

### 化标准型解题步骤（必背）

| 步骤 | 操作 |
|------|------|
| **Step 1** | 识别系数 $a_{11}, a_{12}, a_{22}$ |
| **Step 2** | 计算判别式 $\Delta = a_{12}^2 - a_{11}a_{22}$，判断类型 |
| **Step 3** | 写出特征方程 $a_{11}(dy)^2 - 2a_{12}dxdy + a_{22}(dx)^2 = 0$ |
| **Step 4** | 求解特征方程，得到特征曲线 $\varphi(x,y)=c_1, \psi(x,y)=c_2$ |
| **Step 5** | 作变量替换 $\xi = \varphi(x,y), \eta = \psi(x,y)$ |
| **Step 6** | 利用复合函数求导法则，计算 $u_x, u_y, u_{xx}, u_{xy}, u_{yy}$ |
| **Step 7** | 代入原方程，化简得到标准型 |

---

### 真题1（2023级期末试题）⭐⭐⭐⭐⭐

**题目**：化简二阶线性偏微分方程为标准型

$$
x^2u_{xx} + 2xyu_{xy} + y^2u_{yy} + xyu_x + y^2u_y = 0
$$

**解答：**

**Step 1. 识别系数**
$a_{11} = x^2, \quad a_{12} = xy, \quad a_{22} = y^2$

**Step 2. 判断类型**
$$
\Delta = a_{12}^2 - a_{11}a_{22} = (xy)^2 - x^2 \cdot y^2 = 0
$$
故方程为**抛物型**。

**Step 3. 特征方程**
$$
x^2(dy)^2 - 2xydxdy + y^2(dx)^2 = 0
$$
即
$$
(xdy - ydx)^2 = 0 \implies xdy = ydx
$$

**Step 4. 求解特征曲线**
$$
\frac{dy}{y} = \frac{dx}{x} \implies \ln|y| = \ln|x| + C \implies \frac{y}{x} = C
$$
令 $\xi = \frac{y}{x}$，再取 $\eta = x$（使得雅可比行列式不为零）。

**Step 5. 计算偏导数**
$$
\xi_x = -\frac{y}{x^2}, \quad \xi_y = \frac{1}{x}
$$
$$
\eta_x = 1, \quad \eta_y = 0
$$

一阶偏导：
$$
u_x = u_\xi \xi_x + u_\eta \eta_x = -\frac{y}{x^2}u_\xi + u_\eta
$$
$$
u_y = u_\xi \xi_y + u_\eta \eta_y = \frac{1}{x}u_\xi
$$

二阶偏导：
$$
u_{xx} = \frac{\partial}{\partial x}\left(-\frac{y}{x^2}u_\xi + u_\eta\right) = \frac{2y}{x^3}u_\xi + \frac{y^2}{x^4}u_{\xi\xi} - \frac{2y}{x^2}u_{\xi\eta} + u_{\eta\eta}
$$
$$
u_{xy} = \frac{\partial}{\partial y}\left(-\frac{y}{x^2}u_\xi + u_\eta\right) = -\frac{1}{x^2}u_\xi - \frac{y}{x^3}u_{\xi\xi} + \frac{1}{x}u_{\xi\eta}
$$
$$
u_{yy} = \frac{\partial}{\partial y}\left(\frac{1}{x}u_\xi\right) = \frac{1}{x^2}u_{\xi\xi}
$$

**Step 6. 代入原方程化简**
代入后可得标准型：
$$
u_{\eta\eta} + \frac{1}{\eta}u_\eta = 0
$$

---

### 真题2（2020级期末试题）⭐⭐⭐⭐⭐

**题目**：把方程

$$
u_{xx} - 2\lambda u_{xy} - 3\lambda^2 u_{yy} + u_x + \lambda u_y = 0
$$

化为标准型，并且求出通解，其中 $\lambda$ 为常数。

**解答：**

**Step 1. 识别系数**
$a_{11} = 1, \quad a_{12} = -\lambda, \quad a_{22} = -3\lambda^2$

**Step 2. 判断类型**
$$
\Delta = a_{12}^2 - a_{11}a_{22} = \lambda^2 - 1 \cdot (-3\lambda^2) = 4\lambda^2 \geq 0
$$

**分情况讨论：**

---

**情形1：$\lambda \neq 0$ 时，$\Delta > 0$，双曲型**

特征方程：
$$
(dy)^2 + 2\lambda dxdy - 3\lambda^2(dx)^2 = 0
$$
$$
(dy + 3\lambda dx)(dy - \lambda dx) = 0
$$

解得：
$$
dy = -3\lambda dx \quad \text{或} \quad dy = \lambda dx
$$

积分得特征曲线：
$$
y + 3\lambda x = c_1, \quad y - \lambda x = c_2
$$

令 $\xi = y + 3\lambda x, \quad \eta = y - \lambda x$

计算偏导数代入原方程，化简可得标准型：
$$
u_{\xi\eta} = \frac{1}{4\lambda}u_\xi
$$

进一步求通解：
令 $v = u_\xi$，则 $v_\eta = \frac{1}{4\lambda}v$
解得：$v = f(\xi)e^{\eta/(4\lambda)}$
因此：$u = \int f(\xi)e^{\eta/(4\lambda)}d\xi + g(\eta)$
最终通解：
$$
u(\xi,\eta) = F(\xi)e^{\eta/(4\lambda)} + g(\eta)
$$
或
$$
u(x,y) = F(y + 3\lambda x)e^{(y - \lambda x)/(4\lambda)} + g(y - \lambda x)
$$

---

**情形2：$\lambda = 0$ 时，$\Delta = 0$，抛物型**

原方程简化为：
$$
u_{xx} + u_x = 0
$$

这已经是标准型。

求通解：令 $v = u_x$，则 $v_x + v = 0$
解得：$v = f(y)e^{-x}$
因此：
$$
u(x,y) = \int f(y)e^{-x}dx + g(y) = -f(y)e^{-x} + g(y)
$$

---

## 本章小结

| 类型 | 判别式条件 | 标准形式 | 典型方程 | 考试频率 |
|------|------------|----------|----------|----------|
| **双曲型** | $\Delta > 0$ | $u_{ss} - u_{tt} + \cdots = 0$ | 波动方程 | ⭐⭐⭐ |
| **抛物型** | $\Delta = 0$ | $u_{\eta\eta} + \cdots = 0$ | 热传导方程 | ⭐⭐⭐⭐⭐ |
| **椭圆型** | $\Delta < 0$ | $u_{\alpha\alpha} + u_{\beta\beta} + \cdots = 0$ | Laplace方程 | ⭐⭐ |

> 💡 **考试提醒**：化标准型题目每年必考12分！务必掌握解题步骤，特别注意抛物型的化简计算！

---

> 本文件基于PPT内容整理优化，补充了历年考试真题和解题模板。

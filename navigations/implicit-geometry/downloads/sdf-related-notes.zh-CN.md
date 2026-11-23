# 《Variational Analysis on the Signed Distance Functions》中文解读

作者：Honglin Luo、Xianfu Wang、Brett Lukens。发表于 *Journal of Optimization Theory and Applications* 180 (2019), 751–774；2018 年在线发表。DOI：10.1007/s10957-018-1414-2。

来源为同目录 `sdf-related.pdf`，共 24 页。下文“PDF 第 p 页”对应期刊页码 750+p。保留原文结果编号；教学例子与额外推导会明确标出。本文的 SDF 指 Signed Distance Function（有符号距离函数）。

## 1. 文章研究什么

文章研究如何用一个实值函数表达集合的内外关系、边界距离和法向几何，并将这些几何信息与凸分析、次微分、单调算子联系起来。

其主线为：

1. 一般集合：SDF 的定义、集合辨识、几何变换，以及 Eikonal 方程几乎处处成立的精确条件。
2. 闭凸集合：最近边界点、法向分解、骨架，以及边界投影的参数单调性与极大单调性。
3. 凸分析：SDF 是支撑超平面有符号距离的上确界，因而是凸函数；由此研究 Fenchel 共轭和次微分。
4. 构造应用：利用 SDF 构造“函数是凸的，但次微分定义域不是凸集”的例子。
5. 显式计算：给出非负正交锥、二阶锥的 SDF 和边界投影。

文章主要贡献是理论分析。碰撞检测、分类、形状分析、水平集等是其动机；文章没有为这些任务给出新的实验算法。

## 2. 符号与最基本定义

设所有距离均使用有限维欧氏空间的二范数。通常只考虑非空且不等于整个空间的集合，称为 proper set。

| 符号 | 含义 |
|---|---|
| \(S^c\) | \(\mathbb R^n\setminus S\)，补集 |
| \(\operatorname{cl}S,\operatorname{int}S,\partial S\) | 闭包、内部、边界；此处 \(\partial S\) 等同原文 bdry S |
| \(\mathbb B,\mathbb S\) | 闭单位球、单位球面 |
| \(\operatorname{conv}A,\operatorname{ri}A\) | 凸包、相对内部 |
| \(d_S(x)=\inf_{s\in S}\|x-s\|\) | 普通距离函数 |
| \(P_S(x)=\{s\in S:\|x-s\|=d_S(x)\}\) | 到集合的投影，可能多值或为空 |
| \(Q_S(x)=P_{\partial S}(x)\) | 到边界的投影；proper 集合的边界闭且非空，因此始终非空 |
| \(\sigma_C(v)=\sup_{c\in C}\langle v,c\rangle\) | 支撑函数，允许值为 \(+\infty\) |
| \(\iota_A\) | 集合示性函数：在 A 上为 0，在外部为 \(+\infty\) |
| \(N_C(c)\) | 外法锥：\(\{v:\langle v,y-c\rangle\le0,\ \forall y\in C\}\)，当 \(c\in C\) |
| \(f^*(v)=\sup_x\{\langle v,x\rangle-f(x)\}\) | Fenchel 共轭 |
| \(\partial f(x)\) | 凸次微分：所有满足 \(f(y)\ge f(x)+\langle v,y-x\rangle\) 的 v |
| \(\operatorname{dom}\partial f\) | 有至少一个有限次梯度的点集 |
| \(\operatorname{ran}\partial f\) | 所有点处次梯度的并集 |
| \(K^{\ominus}\) | 负极锥：\(\{v:\langle v,k\rangle\le0,\ \forall k\in K\}\) |

### 2.1 SDF 定义：式 (1)，PDF 第 3 页

\[
D_S(x)=\begin{cases}
d_S(x),&x\notin S,\\
-d_{S^c}(x),&x\in S.
\end{cases}
\]

即外正、内负、边界为零。负数的绝对值表示内部点离开集合至少要移动多远。

Fact 1（PDF 第 4 页）给出：

\[
D_S=d_S-d_{S^c},\qquad
|D_S|=d_S+d_{S^c}=d_{\partial S},\qquad
D_S^2=d_S^2+d_{S^c}^2.
\]

\[
\operatorname{int}S=\{D_S<0\},\quad
\operatorname{cl}S=\{D_S\le0\},\quad
\partial S=\{D_S=0\},\quad
\operatorname{int}S^c=\{D_S>0\}.
\]

注意：一般集合满足的是 \(\{D_S\le0\}=\operatorname{cl}S\)。只有集合闭时，才可直接写成 \(S=\{D_S\le0\}\)。

**例子 A：半径为 2 的闭圆盘（教学例子，原文单位球例子的缩放）。**

\[
C=\{x:\|x\|\le2\},\quad D_C(x)=\|x\|-2,\quad d_C(x)=\max\{\|x\|-2,0\}.
\]

| 点 | SDF | 普通距离 | 意义 |
|---|---:|---:|---|
| (3,0) | 1 | 1 | 外部距边界 1 |
| (1,0) | -1 | 0 | 内部深度 1 |
| (0,0) | -2 | 0 | 中心深度 2 |
| (2,0) | 0 | 0 | 边界 |

普通距离在整个集合内部都是 0；SDF 保留了内部深度。

### 2.2 一般集合的基本运算：命题 3.1–3.4

\[
D_{S^c}=-D_S,\qquad S\subseteq R\Longrightarrow D_R\le D_S.
\]

\[
D_{S+w}(x)=D_S(x-w),\qquad D_{kS}(x)=kD_S(x/k)\quad(k>0).
\]

若 \(A A^T=k^2I\)，则 A 是旋转或反射加均匀缩放，且

\[
D_{AS}(x)=kD_S(A^{-1}x)=kD_S(k^{-2}A^Tx).
\]

若 \(\phi\) 为欧氏等距变换，则 \(D_{\phi(S)}=D_S\circ\phi^{-1}\)。这些恒等式不直接适用于任意非均匀缩放。

若 \(\operatorname{int}S=\varnothing\)，则 \(D_S=d_S\)。对闭的 proper 集合，若 \(D_R=\alpha D_S\) 且 \(\alpha>0\)，则 \(R=S,\alpha=1\)：同一个零水平集上的任意缩放通常已不是精确距离函数。

对闭的 proper 集合：\(D_S(\alpha x)=\alpha D_S(x)\) 当且仅当 \(\alpha S=S\)。若 S 还凸，则对所有正比例均齐次，当且仅当 S 是闭凸锥。

**应用：** 刚体运动后直接反变换查询点；均匀缩放模型时同时缩放距离值；检查一个隐式函数是否仍是精确欧氏 SDF。

### 2.3 SDF 到底区分哪些集合：命题 3.5，PDF 第 6 页

\[
D_R=D_S
\iff\operatorname{cl}R=\operatorname{cl}S\text{ 且 }\operatorname{int}R=\operatorname{int}S.
\]

也等价于闭包相同且补集闭包相同，或闭包相同且边界相同。这是 Brezis–Haraux 的“几乎相等”（拓扑意义），不等于“只差一个零测集”。

**例子 B：开球与闭球。** 它们的闭包、内部相同，故 SDF 完全相同；SDF 不记录边界点是否属于集合。

**例子 C：原文例 3.1。** 令 S 为单位圆盘中的有理坐标点，则

\[
D_S(x)=\max\{\|x\|-1,0\},\quad
D_{\operatorname{cl}S}(x)=\|x\|-1.
\]

S 虽稠密却无内部。普通距离看不出它与闭圆盘的差别，SDF 能看出。也说明一般非闭集合可能不凸，但 SDF 却凸。

推论 3.1 将此联系到算子值域：对 proper、下半连续凸函数 f,g，若 \(\operatorname{dom}f\cap\operatorname{int}\operatorname{dom}g\ne\varnothing\)，则

\[
D_{\operatorname{ran}(\partial f+\partial g)}
=D_{\operatorname{ran}\partial f+\operatorname{ran}\partial g}.
\]

这里右端下标的加法是 Minkowski 和。结论是两个值域集合具有相同 SDF，并非直接断言它们逐点相等。

### 2.4 Lipschitz 性与 Eikonal 方程：定理 3.1

对一般 proper 集合，SDF 是 1-Lipschitz：

\[
|D_S(x)-D_S(y)|\le\|x-y\|.
\]

这是背景性质，本文在证明中使用；补充解释如下：同号时使用距离函数的 Lipschitz 性；异号时线段必穿过边界，两个端点到边界的距离之和不超过线段长度。因此它几乎处处可微，但仅此不能推出梯度模为 1。

**定理 3.1（PDF 第 6 页）：**

\[
\left[D_S\text{ 几乎处处可微且 }\|\nabla D_S\|=1\text{ 几乎处处}\right]
\iff\mu(\partial S)=0.
\]

\(\mu\) 是 n 维 Lebesgue 测度。对 proper 凸集，该条件自动成立（推论 3.2）。

在非边界的可微点，若最近边界点为 q，则

\[
\nabla D_S(x)=\frac{x-q}{D_S(x)}.
\]

内部 D 为负，因此梯度仍指向外部。梯度不存在的点不能套用这个单值公式。

**例子 D：** 圆盘中心不可微，其余点梯度为单位径向向量；一个点集的不可微性不妨碍“几乎处处”结论。相反，取 \(S=\mathbb Q^n\)，其边界为整个空间，\(D_S\equiv0\)，梯度处处为 0，不能满足单位梯度条件。这是教学反例。

**应用：** 水平集中的距离校准、距离函数近似的梯度检查。若用 \((\|\nabla f\|-1)^2\) 作为学习或数值误差指标，这是由理论启发的应用，本文并未提出学习算法；单位梯度也不要求在中心、棱角、骨架上逐点成立。

## 3. 边界投影、法向分解和骨架

本节除特别说明外，C 均为非空、proper、闭凸集。

### 3.1 普通投影与边界投影：式 (2)

\[
Q_S(x)=P_{\partial S}(x)=
\begin{cases}
P_{\operatorname{cl}S}(x),&x\in S^c,\\
P_{\operatorname{cl}(S^c)}(x),&x\in S.
\end{cases}
\]

对闭凸 C，外部点的 \(Q_C=P_C\) 且单值；内部点的 \(P_C(x)=x\)，而 \(Q_C(x)\) 可能包含多个点。

**例子 E：单位圆盘。** \(x\ne0\) 时 \(Q_C(x)=\{x/\|x\|\}\)；中心 \(Q_C(0)=\mathbb S\)。内部投影到集合本身与投影到边界完全不同。

### 3.2 定理 3.2：沿法线保持同一投影（PDF 第 7 页）

若 \(q\in Q_C(x)\)：

\[
x\in\operatorname{int}C,\ \alpha>0
\Longrightarrow Q_C((1-\alpha)x+\alpha q)=\{q\}.
\]

\[
x\notin\operatorname{int}C,\ \alpha\ge0
\Longrightarrow Q_C(\alpha x+(1-\alpha)q)=\{q\}.
\]

从内部点向选定的最近边界点移动，越过边界后沿该方向继续走，最近点仍是 q。注意第一式排除了 \(\alpha=0\)，因为初始点可能本来就有多个最近点。

命题 3.9 进一步给出沿这些射线的线性变化：

\[
D_C((1-\alpha)x+\alpha q)=(1-\alpha)D_C(x)
\quad(x\in\operatorname{int}C,\alpha\ge0),
\]

\[
D_C(\alpha x+(1-\alpha)q)=\alpha D_C(x)
\quad(x\notin\operatorname{int}C,\alpha\ge0).
\]

### 3.3 定理 3.3：统一的内外法向分解（PDF 第 8 页）

\[
\boxed{x=q+D_C(x)n,\quad q\in Q_C(x),\quad n\in\mathbb S\cap N_C(q).}
\]

SDF 是带符号的法向坐标。命题 3.6 说明内部时 \(q-x\in N_C(q)\)，外部时 \(x-q\in N_C(q)\)。

可微时可写为

\[
q=x-D_C(x)\nabla D_C(x).
\]

这是上述分解与梯度公式的直接推论，不是原文另外编号的定理。

**例子 F：半径 2 的圆盘。** 对 \(x=(1,0)\)，\(q=(2,0),D=-1,n=(1,0)\)，故 \((1,0)=(2,0)-1(1,0)\)。对 \(x=(3,0)\)，同一 q,n 配合 \(D=1\)。

**应用：** 将空间点还原为边界点和法向位移；从 SDF 获得接触点、单位外法向；做边界附近的局部坐标计算。多值点需要保留集合信息。

### 3.4 骨架与定理 3.4

\[
\operatorname{skel}(C)=\{x\in C:Q_C(x)\text{ 不为单点}\}.
\]

命题 3.7：骨架位于内部且没有内部；两个不同骨架点的最近边界点集合不相交。

引理 3.3：

\[
\operatorname{skel}(C)=\varnothing
\iff C\text{ 是半空间，或 }\operatorname{int}C=\varnothing.
\]

**定理 3.4（PDF 第 9 页）：** 若骨架非空、\(x\in\operatorname{int}C,q\in Q_C(x)\)，设 \(x_\gamma=x+\gamma(q-x)\)，则

\[
-\infty<\beta:=\min\{\gamma:q\in Q_C(x_\gamma)\}\le0.
\]

若 x 本身为骨架点，则 \(\beta=0\)。也就是说，沿选定投影方向反向深入时，该投影能保持为最近点的区间有一个有限端点。

**例子 G：单位圆盘。** \(x=(0.4,0),q=(1,0)\)，则 \(x_\gamma=(0.4+0.6\gamma,0)\)，端点 \(\beta=-2/3\) 对应圆心。越过圆心后 q 不再最近。

原文又定义 V(C)：对每一个最近边界点，向反方向移动立即失去该最近点性质的内部点集。引理 3.4 给出每个内部点的唯一 V(C) 起点分解：

\[
x=y+\alpha(q-y),\quad y\in V(C),\ q\in Q_C(y),\quad0\le\alpha<1.
\]

文章只证明 \(\operatorname{skel}(C)\subseteq V(C)\)，并明确说未能证明二者一般相等。不要把 V(C) 与骨架无条件视为同一个集合。

### 3.5 定理 3.5–3.6：单调算子性质

**单调性：** 若 \(p\in P_S(x),q\in P_S(y)\)，则 \(\langle x-y,p-q\rangle\ge0\)。

**参数单调性（paramonotonicity）：** 若上式等号成立，则还可交换投影值：\(q\in P_S(x),p\in P_S(y)\)。

**定理 3.5（PDF 第 10 页）：** 任意非空 proper 闭集 S 的多值投影 \(P_S\) 都参数单调，无须 S 凸。因此 \(Q_C=P_{\partial C}\) 参数单调（推论 3.3）。

证明核心就是最近距离不等式：

\[
\|x-p\|^2+\|y-q\|^2
\le\|x-q\|^2+\|y-p\|^2,
\]

右端减左端等于 \(2\langle x-y,p-q\rangle\)。若等号成立，两个距离极小值都能交换达到。

**极大单调：** 算子图像已不能添加新的输入输出对而保持单调。

**定理 3.6（PDF 第 10–11 页）：**

\[
Q_C\text{ 极大单调}\iff\operatorname{skel}(C)=\varnothing.
\]

**例子 H：单位圆盘。** 中心的 \(Q_C(0)=\mathbb S\) 不是凸集，因而无法成为极大单调算子的该点值；半空间的边界是一张超平面，其边界投影则极大单调。参数单调并不保证单值或连续。

### 3.6 边界候选点：式 (5)、命题 3.8

\[
Z_C(x)=\begin{cases}
\{c\in\partial C:x=c-z,\ z\in N_C(c)\},&x\in\operatorname{int}C,\\
\{c\in\partial C:x=c+z,\ z\in N_C(c)\},&x\notin\operatorname{int}C.
\end{cases}
\]

\[
Q_C(x)=\operatorname*{argmin}_{c\in Z_C(x)}\|x-c\|.
\]

外部或边界的 \(Z_C(x)\) 是单点，内部可能有多个候选，法向条件只是候选条件，仍需比较距离。

**例子 I：单位圆盘内 \(x=(0.5,0)\)。** \((1,0)\) 和 \((-1,0)\) 都满足向内法向候选条件，但距离分别为 0.5 和 1.5，只有前者属于 Q。

## 4. SDF 为什么是凸函数：定理 3.7

### 4.1 支撑超平面的表示（PDF 第 11–12 页）

闭凸集由所有支撑半空间的交表示：

\[
C=\bigcap_{\|v\|=1}H(v),\quad H(v)=\{x:\langle v,x\rangle\le\sigma_C(v)\}.
\]

\(\sigma_C(v)=+\infty\) 的方向可忽略。单位法向半空间的 SDF 是仿射函数：

\[
D_{H(v)}(x)=\langle v,x\rangle-\sigma_C(v).
\]

**定理 3.7：**

\[
\boxed{D_C(x)=\sup_{\|v\|=1}\{\langle v,x\rangle-\sigma_C(v)\}
=(\sigma_C+\iota_{\mathbb S})^*(x).}
\]

凸性的原因非常直接：仿射函数的逐点上确界是凸函数。

内部时，每个式子都是到一张支撑平面的负距离，最大者对应最近的平面；外部时，投影点处的支撑平面提供恰好等于距离的值。注意上确界遍历所有单位方向，而非只遍历多面体已有约束的法向。

### 4.2 用圆盘推导一次

若 \(C=R\mathbb B\)，则 \(\sigma_C(v)=R\|v\|\)。因此

\[
D_C(x)=\sup_{\|v\|=1}\langle v,x\rangle-R=\|x\|-R.
\]

这是 Cauchy–Schwarz 不等式的直接应用。对 x 非零，最优方向为 \(v=x/\|x\|\)。

### 4.3 共轭：推论 3.4、式 (10)

\[
D_C^*=(\sigma_C+\iota_{\mathbb S})^{**}
=\operatorname{cl}\operatorname{conv}(\sigma_C+\iota_{\mathbb S}).
\]

这里对函数取闭凸包，等价于对其上图集取闭凸化；不是简单对数值或对 C 取凸包。

**圆盘例子的直接推导：**

\[
D_{R\mathbb B}^*(v)=R+\iota_{\mathbb B}(v).
\]

如果 \(\|v\|\le1\)，\(\langle v,x\rangle-\|x\|+R\le R\) 且 x=0 达到 R；若 \(\|v\|>1\)，沿 v 方向令 x 趋于无穷就得到 \(+\infty\)。

### 4.4 凸集条件的边界

SDF 的凸性对非闭凸集和 nearly convex 集合也成立（Remark 3.1、推论 3.5）。Nearly convex 表示存在凸集 A 满足 \(A\subseteq C\subseteq\operatorname{cl}A\)，此时 \(D_C=D_{\operatorname{cl}C}=D_{\operatorname{ri}C}\)。

一般集合的 SDF 凸，不必推出该集合凸，原文有理点圆盘即是反例。但若已知 C 是闭集，则反向成立，因为 \(C=\{D_C\le0\}\) 是凸函数的次水平集。这是直接推论。

### 4.5 多面体：推论 3.6、式 (11)

设 \(C=\bigcap_{i=1}^m\{x:a_i^Tx\le b_i\}\) 非空，且 \(\|a_i\|=1\)，则

\[
D_C(x)=\begin{cases}
\max_i(a_i^Tx-b_i),&x\in C,\\
d_C(x),&x\notin C.
\end{cases}
\]

外部距离可先解二次规划

\[
y^*\in\operatorname*{argmin}_{a_i^Ty\le b_i}\frac12\|x-y\|^2,
\quad d_C(x)=\|x-y^*\|.
\]

**例子 J：正方形 \(C=[-1,1]^2\)。** 内部 \(D_C(x,y)=\max(|x|,|y|)-1\)；例如 \(D_C(0.2,0.7)=-0.3\)。但外部 \((2,2)\) 的真实距离为 \(\sqrt2\)，仅取约束最大违反量得到 1，并不准确。

全空间的精确式为（教学推导）

\[
D_C(x,y)=\sqrt{\max(|x|-1,0)^2+\max(|y|-1,0)^2}
+\min\{\max(|x|-1,|y|-1),0\}.
\]

**应用：** 用 \(D_C\le0\) 编码闭凸可行域；用负值衡量可行点的安全余量；外部最近可行点通过凸二次规划计算。上面的有限 max 公式只能在内部当作精确 SDF。

## 5. 次微分与最大内接球：定理 3.8–3.9

### 5.1 次微分的几何意义

在光滑处次微分只包含梯度；在不可微处，它收集所有支撑斜率。对于本节凸 SDF，Clarke 次微分与 Fenchel 次微分一致。

**定理 3.8（PDF 第 13–14 页）：**

\[
\boxed{\partial D_C(x)=\begin{cases}
\left\{\dfrac{x-P_C(x)}{D_C(x)}\right\},&x\notin C,\\[6pt]
\operatorname{conv}(\mathbb S\cap N_C(x)),&x\in\partial C,\\[6pt]
\left\{\dfrac{x-q}{D_C(x)}:q\in\operatorname{conv}Q_C(x)\right\},&x\in\operatorname{int}C.
\end{cases}}
\]

第三行也等于 \(\operatorname{conv}\{(x-q)/D_C(x):q\in Q_C(x)\}\)。内部有多个最近边界点，便产生多组候选外法向，再取凸包。

**单位球例子（原文例 3.4、3.6）：**

\[
\partial D_{\mathbb B}(x)=\begin{cases}
\{x/\|x\|\},&x\ne0,\\
\mathbb B,&x=0.
\end{cases}
\]

圆心不是梯度为零的光滑点；它没有梯度，但次微分包含 0。

**正方形中心（教学例子）：**

\[
Q_C(0)=\{(1,0),(-1,0),(0,1),(0,-1)\},
\quad\partial D_C(0)=\{v:|v_1|+|v_2|\le1\}.
\]

这里是菱形，不是欧氏单位圆盘；局部次梯度集合会反映形状。

**切勿混淆普通距离与 SDF 的边界次微分。** 普通距离为

\[
\partial d_C(x)=\mathbb B\cap N_C(x)\quad(x\in C),
\]

而 SDF 边界为 \(\operatorname{conv}(\mathbb S\cap N_C(x))\)。例如单位球边界点 q，\(\partial D_C(q)=\{q\}\)，但 \(\partial d_C(q)=[0,q]\)。当 C 无内部时二者才整体相同（原文例 3.3）。

### 5.2 定理 3.9：SDF 的极小点（PDF 第 15 页）

若 C 无内部，则 \(\operatorname*{argmin}D_C=C\)。若 C 有内部，则

\[
\boxed{x\in\operatorname*{argmin}D_C
\iff x\in\operatorname{int}C\text{ 且 }x\in\operatorname{conv}Q_C(x).}
\]

由凸最优性条件 \(0\in\partial D_C(x)\) 直接得到。这些点是最大内接球的中心；最小值的负数为内接球最大半径。推论 3.7 说明这些极小点必在骨架上。

**例子 K：长方形 \([-2,2]\times[-1,1]\)（教学例子）。** 最大内接圆半径为 1，所有 \((s,0),|s|\le1\) 都是圆心。极小点无需唯一。上下最近边界点的中点就是圆心，满足凸包条件。

半空间的 SDF 是仿射函数，向内部走可无限减小，因而没有极小点；闭凸且有内部本身并不保证最大内接球存在。

**应用：** 求安全中心、几何裕量最大的可行点、最大内接球；用次梯度方法处理中心和骨架的不可微性。

## 6. 次微分值域：定理 3.10–3.13

定义

\[
A=\mathbb S\cap\operatorname{dom}\sigma_C,
\qquad B=\mathbb S\cap\operatorname{ran}N_C.
\]

这里 B 是本节临时记号，不是单位球 \(\mathbb B\)。dom 支撑函数表示支撑值有限；ran 法锥表示某处确实存在该支撑法向。无界集合可能有有限但未达到的支撑值，故这两类集合不宜随意混同。

### 6.1 定理 3.10（PDF 第 15 页）

\[
\operatorname{ri}\operatorname{conv}A
\subseteq\operatorname{ran}\partial D_C
\subseteq\operatorname{cl}\operatorname{conv}A
\subseteq\mathbb B,
\]

\[
\operatorname{ri}\operatorname{ran}\partial D_C
=\operatorname{ri}\operatorname{conv}A.
\]

### 6.2 定理 3.11（PDF 第 16 页）

\[
\operatorname{ri}\operatorname{conv}B
\subseteq\operatorname{ran}\partial D_C
\subseteq\operatorname{conv}B
\subseteq\mathbb B,
\]

\[
\operatorname{ri}\operatorname{ran}\partial D_C
=\operatorname{ri}\operatorname{conv}B.
\]

所有次梯度的模都不超过 1，与 1-Lipschitz 性一致。值域是 nearly convex，但未必闭或凸；边界点是否存在需要单独判断。

推论 3.8：C 若紧且凸，则 \(\operatorname{ran}\partial D_C=\mathbb B\)。单位球中心一个点就贡献整个单位球；其他形状可能需要合并多个点处的次梯度。

推论 3.9：对 proper 闭凸锥 K，

\[
\operatorname{ran}\partial D_K=\operatorname{conv}(\mathbb S\cap K^{\ominus}).
\]

### 6.3 定理 3.12：普通距离的值域（PDF 第 16 页）

\[
\operatorname{ran}\partial d_C
=\partial d_C(\partial C)
=\bigcup_{x\in\partial C}(\mathbb B\cap N_C(x))
=\mathbb B\cap\operatorname{ran}N_C.
\]

若 C 紧，则等于 \(\mathbb B\)；若 C=K 为闭凸锥，则等于 \(\mathbb B\cap K^{\ominus}\)。C 无内部时，SDF 与普通距离相同，值域也相同（推论 3.10）。

**半空间例子：** 单位外法向为 a，则 \(\operatorname{ran}\partial D_H=\{a\}\)，而 \(\operatorname{ran}\partial d_H=[0,a]\)。这直接说明两种值域公式不能互换。

### 6.4 定理 3.13（PDF 第 16–17 页）

\[
\operatorname{ran}\partial D_C
=\partial D_C(\partial C\cup\operatorname{skel}(C)).
\]

右侧表示在这些点处取次微分的并集。若只研究所有可能的次梯度，观察边界与骨架就够了。非骨架内部点和外部点贡献的单个单位法向，已在边界次微分中出现。

推论 3.11：

\[
\partial D_C(C^c)=\mathbb S\cap\operatorname{ran}N_C
\subseteq\operatorname{ran}\partial D_C.
\]

原文例 3.5：\(C=\operatorname{epi}(x\mapsto x^2)\) 满足

\[
N_C(x,x^2)=\mathbb R_+(2x,-1),\quad
\operatorname{ran}N_C=\{(u,v):v<0\}\cup\{(0,0)\}.
\]

法向可以无限接近水平，却不存在非零水平支撑法向，说明法锥值域不必闭，相关公式中的闭包不可随意删除。

## 7. 用 SDF 构造非凸的次微分定义域

### 7.1 为什么“凸函数”可以有“非凸次微分定义域”

凸函数的有效定义域 \(\operatorname{dom}f\) 是凸集，但边界点处可能没有有限的支撑斜率，因此 \(\operatorname{dom}\partial f\) 可以缺少部分边界点。两个概念需要区分。

设 \(\Gamma_0\) 是 proper、下半连续凸函数类，\(\Lambda_n=\{\operatorname{dom}\partial f:f\in\Gamma_0(\mathbb R^n)\}\)。

### 7.2 命题 4.3、定理 4.1（PDF 第 18–19 页）

给定有限凸函数 g，令

\[
F(t,x)=\begin{cases}
\max\{g(x),-\sqrt t\},&t\ge0,\\
+\infty,&t<0.
\end{cases}
\]

则 F 是 proper、下半连续凸函数。若 \(M=\{x:g(x)<0\}\ne\varnothing\)，则

\[
\operatorname{dom}F=\mathbb R_+\times\mathbb R^n,
\quad\operatorname{dom}\partial F
=(\mathbb R_+\times\mathbb R^n)\setminus(\{0\}\times M).
\]

选 \(g=D_C\)，其中 C 为有内部的 proper 闭凸集，就得到定理 4.1：

\[
\boxed{(\mathbb R_+\times\mathbb R^n)
\setminus(\{0\}\times\operatorname{int}C)\in\Lambda_{n+1}.}
\]

**为什么边界这部分没有次梯度？** 当 \(D_C(x)<0\) 时，固定 x 并令 t 从 0 向右移动一小段，F 由 \(-\sqrt t\) 主导。其右导数趋于 \(-\infty\)，不存在有限支撑斜率。具体地，若在 (0,x) 存在 t 方向次梯度 a，则须 \(-\sqrt t\ge at\) 对所有足够小的正 t 成立，即 \(a\le-1/\sqrt t\)，任何有限 a 都做不到。

### 7.3 具体的非凸例子（教学展开）

取 \(C=[-1,1]\)，\(D_C(x)=|x|-1\)，则

\[
F(t,x)=\max\{|x|-1,-\sqrt t\}\quad(t\ge0),
\]

\[
\operatorname{dom}\partial F
=\{(t,x):t>0\}\cup\{(0,x):|x|\ge1\}.
\]

\((0,-1)\) 与 \((0,1)\) 在该集合中，但中点 \((0,0)\) 不在，故次微分定义域不是凸集。其内部仍是 \(t>0\)，闭包仍是 \(t\ge0\)，因而 nearly convex。

### 7.4 定理 4.2：从已有定义域上删除一块边界

若 \(E\in\Lambda_{n+1}\)，半空间 H 包含 E，且 \(\operatorname{ri}H\cap\operatorname{ri}E\ne\varnothing\)；若 \(C\subsetneq\partial H\) 为闭凸集且 \(\dim\operatorname{aff}C=n\)，则

\[
E\setminus\operatorname{ri}C\in\Lambda_{n+1}.
\]

证明使用命题 4.1：有限个次微分定义域在相对内部具有公共点时，其交仍可作为次微分定义域。推论 4.1–4.4 与命题 4.2 允许做可逆线性变换、平移、等距变换，使删除操作适用于不同朝向的半空间。

**用途：** 构造凸分析与单调算子理论中的反例；理解优化模型为何“函数有定义”却不一定“存在有限次梯度”。这是文章的理论应用，不是碰撞检测算法。

## 8. 凸锥的显式公式：定理 5.1 和例 5.1–5.2

### 8.1 定理 5.1：闭凸锥的 SDF 是次线性的（PDF 第 20 页）

对 proper 闭凸锥 K，令 \(A=\operatorname{conv}(K^{\ominus}\cap\mathbb S)\)，则

\[
D_K=(\iota_{K^{\ominus}\cap\mathbb S})^*
=\sigma_{K^{\ominus}\cap\mathbb S}=\sigma_A,
\quad D_K^*=\iota_A.
\]

K 的支撑函数为 \(\sigma_K=\iota_{K^{\ominus}}\)，将其代入定理 3.7 即得。

\[
D_K(\alpha x)=\alpha D_K(x)\ (\alpha\ge0),\quad
D_K(x+y)\le D_K(x)+D_K(y),\quad |D_K(x)|\le\|x\|.
\]

“次线性”即正齐次且次可加，不表示函数只能为非负。

命题 5.1 的一般版本为

\[
D_C(x+y)\le\inf_{0<\lambda<1}
\{D_{(1-\lambda)C}(x)+D_{\lambda C}(y)\}.
\]

### 8.2 非负正交锥：原文例 5.1

\[
K=\mathbb R_+^n,\quad
D_K(x)=\begin{cases}
\|(\min(x_1,0),\ldots,\min(x_n,0))\|,&x\notin K,\\
-\min_i x_i,&x\in K.
\end{cases}
\]

若外部，\(Q_K(x)=\{x_+\}\)，逐坐标将负值截断为 0。若内部或边界，令 \(I(x)=\{i:x_i=\min_jx_j\}\)，则

\[
Q_K(x)=\{x-x_i e_i:i\in I(x)\}.
\]

二维数值例子：

| 点 | SDF | 最近边界点 |
|---|---:|---|
| (3,2) | -2 | (3,0) |
| (-3,4) | 3 | (0,4) |
| (-3,-4) | 5 | (0,0) |
| (2,2) | -2 | (0,2)、(2,0) |

最后一个点处不可微，\(\partial D_K(2,2)=\operatorname{conv}\{(-1,0),(0,-1)\}\)。

**应用：** 非负约束的违反程度、内部裕量及投影。多个负坐标必须以欧氏范数组合，不能只取最负坐标的绝对值。

### 8.3 一般底集生成锥：命题 5.2–5.3

设 C 为非空紧凸集，\(K=\operatorname{cone}(C\times\{1\})\)。对 \(y\in C,s>0\)：

\[
N_K(sy,s)=\{(w,-\langle w,y\rangle):w\in N_C(y)\},
\]

\[
N_K(0)=K^{\ominus}
=\{(w,u):u\le-\sigma_C(w)\}.
\]

命题 5.3 将式 (5) 的候选条件具体化。对查询点 (z,t)，非零候选 (sy,s) 需 \(y\in\partial C,w\in N_C(y),s>0\)：

| 查询点的位置 | 候选方程 | 原点是否为候选 |
|---|---|---|
| \((z,t)\in\operatorname{int}K\) | \(s=t-\langle w,y\rangle,\ y=(z+w)/s\) | 当且仅当 \(\sigma_C(-z)\le t\) |
| \((z,t)\notin\operatorname{int}K\cup K^{\ominus}\) | \(s=t+\langle w,y\rangle,\ y=(z-w)/s\) | 否 |
| \((z,t)\in K^{\ominus}\) | 无须其他候选 | 唯一候选为原点 |

内部有多个候选时，仍需按命题 3.8 比较欧氏距离。此工具用于由底集几何推导锥几何。

### 8.4 二阶锥：原文例 5.2，以下补全原文遗漏的分段

\[
K=\{(z,t)\in\mathbb R^n\times\mathbb R:t\ge\|z\|\},\quad r=\|z\|.
\]

\[
\boxed{D_K(z,t)=\begin{cases}
\sqrt{r^2+t^2},&t\le-r,\\
(r-t)/\sqrt2,&t>-r.
\end{cases}}
\]

\[
Q_K(z,t)=\begin{cases}
\{0\},&t\le-r,\\
\{\tfrac12(r+t)(z/r,1)\},&r>0,\ t>-r,\\
\{\tfrac t2(u,1):\|u\|=1\},&r=0,\ t>0.
\end{cases}
\]

为了避免负轴和原点的除零问题，次微分的完整分段可写成：

\[
\partial D_K(z,t)=\begin{cases}
\{(z,t)/\sqrt{r^2+t^2}\},&t\le-r,\ (z,t)\ne0,\\
\{(z/r,-1)/\sqrt2\},&r>0,\ t>-r,\\
\mathbb B_n\times\{-1\}\,/\sqrt2,&r=0,\ t>0,\\
\operatorname{conv}(\mathbb S_{n+1}\cap K^{\ominus}),&(z,t)=0.
\end{cases}
\]

原点处的集合亦为 \(\{(u,v):\|(u,v)\|\le1,\ v\le-1/\sqrt2\}\)。

**数值例子：** 当 \(r=3,t=5\)，点在锥内，\(D=-\sqrt2\)，最近侧面点为 \(4(z/3,1)\)；当 \(r=3,t=1\)，\(D=\sqrt2\)，最近点为 \(2(z/3,1)\)；当 \(r=3,t=-4\)，最近点为原点，\(D=5\)。轴上 \(z=0,t>0\) 具有整圈最近边界点。

**原理：** 利用旋转对称性，将投影化为半径 r 与高度 t 构成的二维截面中，到直线 t=r 或顶点的距离。

**应用：** 二阶锥规划中的可行性距离与投影。上述公式计算的是点到锥边界的带符号距离；锥内普通集合投影仍是点本身。

### 8.5 子空间：原文例 3.8

若 U 为 proper 线性子空间，\(u_1,\ldots,u_k\) 为正交归一基，则

\[
P_U(x)=\sum_i\langle u_i,x\rangle u_i,
\quad Q_U=P_U,\quad D_U=d_U=\|x-P_U(x)\|,
\]

\[
\partial d_U(x)=\begin{cases}
\{(x-P_U(x))/\|x-P_U(x)\|\},&x\notin U,\\
\mathbb B\cap U^\perp,&x\in U,
\end{cases}
\quad\operatorname{ran}\partial d_U=\mathbb B\cap U^\perp.
\]

例如二维中的 x 轴没有二维内部，其 SDF 为 \(|y|\)，不会在轴的一侧取负值。SDF 的“内部深度”针对的是环境空间的内部。

## 9. 原文全部 11 个编号公式索引

| 式号 | PDF 页 | 公式或等价整理 | 作用 |
|---|---:|---|---|
| (1) | 3 | \(D_S=d_S\) 于外部，\(D_S=-d_{S^c}\) 于内部 | 定义 SDF |
| (2) | 4 | \(Q_S=P_{\partial S}\)，按内外分别投影到闭包或补集闭包 | 定义边界投影 |
| (3) | 9 | \(d_{\partial C}(x_\gamma)=\|x_\gamma-q\|=(1-\gamma)\|x-q\|>0\) | 在定理 3.4 反向射线证明中，当 \(\gamma\le0\) 且 q 保持最近点时使用，并非全空间恒等式 |
| (4) | 10 | 距离极小性推出 \(2\langle x-y,p-q\rangle\ge0\)；详细等价式见第 3.5 节 | 证明投影参数单调 |
| (5) | 11 | \(Z_C(x)=\{c\in\partial C:x=c\mp z,z\in N_C(c)\}\)，内部减、外部加 | 生成边界候选点 |
| (6) | 12 | \(D_C(x)=\sup_{\|v\|=1}D_{H(v)}(x)\)，\(x\in C\) | 支撑半空间表示的内部部分 |
| (7) | 12 | \(D_C(x)\ge\sup_{\|v\|=1}D_{H(v)}(x)\)，\(x\notin C\) | 外部一侧不等式 |
| (8) | 12 | \(\sup_{\|v\|=1}D_{H(v)}(x)\ge D_{H(w)}(x)\ge D_C(x)\) | 用分离定理构造方向 w，取得反向不等式 |
| (9) | 12 | \(D_C(x)=\sup_{\|v\|=1}D_{H(v)}(x)\)，\(x\notin C\) | 合并 (7)(8)，完成定理 3.7 |
| (10) | 12 | \(D_C^*=(\sigma_C+\iota_{\mathbb S})^{**}=\operatorname{cl}\operatorname{conv}(\sigma_C+\iota_{\mathbb S})\) | SDF 的 Fenchel 共轭 |
| (11) | 13 | \(\min_{a_i^Ty\le b_i}\frac12\|x-y\|^2\) | 多面体外部投影和距离计算 |

(6)–(9) 是同一个核心表示公式的证明步骤，不是四个互相独立的 SDF 模型。公式中的 \(D^{-1}_C(x)\) 在原文次微分表达式里指 \(1/D_C(x)\)，不是反函数，且只能在 D 非零的分支使用。

## 10. 全部 16 个 Theorem 速查

共同约定：第 3 章涉及 C 的定理除 3.1 与 3.5 外，均假设 C 非空、proper、闭凸。下表列额外条件和核心结论；完整公式见对应正文。

| 定理 | PDF 页 | 核心内容 | 直观例子 / 用途 |
|---|---:|---|---|
| 3.1 | 6 | Eikonal 单位梯度几乎处处成立，当且仅当边界 n 维测度为 0 | 圆盘 vs. 稠密有理点集 |
| 3.2 | 7 | 沿选定最近边界点的法向射线，投影保持同一点并变为单值 | 圆盘径向移动；接触追踪 |
| 3.3 | 8 | \(x=q+D_C(x)n\)，n 是单位外法向 | 用距离与法向恢复边界点 |
| 3.4 | 9 | 骨架非空时，反向法线保持同一投影有有限端点；骨架点端点参数为 0 | 圆盘向内追到圆心 |
| 3.5 | 10 | 任意非空 proper 闭集的多值投影参数单调，不要求凸 | 最近点平方距离比较 |
| 3.6 | 10–11 | 边界投影极大单调，当且仅当骨架为空 | 半空间与圆盘对比 |
| 3.7 | 12 | \(D_C=(\sigma_C+\iota_{\mathbb S})^*\)，因此凸 | 支撑超平面上确界 |
| 3.8 | 13–14 | 外部、边界、内部分别用投影、法锥、最近点凸包计算次微分 | 单位球、正方形中心 |
| 3.9 | 15 | 有内部时，极小点恰为 \(x\in\operatorname{int}C\cap\operatorname{conv}Q_C(x)\) | 最大内接球、长方形中心线段 |
| 3.10 | 15 | 由 \(\mathbb S\cap\operatorname{dom}\sigma_C\) 的凸包夹住次微分值域，并确定相对内部 | 控制可能的次梯度方向 |
| 3.11 | 16 | 改由 \(\mathbb S\cap\operatorname{ran}N_C\) 的凸包夹住值域 | 实际出现的支撑法向 |
| 3.12 | 16 | 普通距离次微分值域为 \(\mathbb B\cap\operatorname{ran}N_C\)，全部可在边界获得 | 与 SDF 值域比较 |
| 3.13 | 16–17 | SDF 次微分的整个值域已在边界与骨架上实现 | 减少几何考察范围 |
| 4.1 | 19 | 从半空间边界删掉 \(\operatorname{int}C\)，仍可成为某凸函数次微分定义域 | \(F=\max(D_C,-\sqrt t)\) |
| 4.2 | 19 | 在相对内部相交等条件下，可从已有次微分定义域删除半空间边界上的凸片相对内部 | 更一般的定义域反例 |
| 5.1 | 20 | 锥 SDF 是极锥单位方向凸包的支撑函数，因而次线性；共轭是该凸包示性函数 | 正交锥、二阶锥 |

## 11. 命题、引理、推论逐条索引

为避免将辅助结果误当作遗漏，下面列出全文其余全部带编号结果。结论中的条件以原文和正文的完整表达为准。

### 命题（Proposition）

| 编号 | PDF 页 | 内容 |
|---|---:|---|
| 3.1 | 4–5 | 补集变号、适当凸性下的闭包/内部不变性、空内部退化、集合辨识、包含关系、空集及全空间的值 |
| 3.2 | 5 | 平移、相似变换、均匀缩放、等距变换的 SDF 变换法则 |
| 3.3 | 5 | 两闭集的 SDF 正比例相等，只能对应同一集合；proper 时比例为 1 |
| 3.4 | 5 | SDF 的缩放齐次性与集合缩放不变性等价；凸情形对应闭凸锥 |
| 3.5 | 6 | SDF 相等的四种等价条件：相同闭包并相同内部/补集闭包/边界 |
| 3.6 | 8 | 内外最近边界点处法锥方向的符号 |
| 3.7 | 8–9 | 骨架无内部；不同骨架点的边界投影集合不相交 |
| 3.8 | 11 | 在法向候选 Z 中选最近者得到 Q；非内部候选唯一 |
| 3.9 | 11 | 法向射线上的 SDF 线性变化 |
| 4.1 | 18 | 有限个次微分定义域在相对内部公共点条件下，对交封闭，且交 nearly convex |
| 4.2 | 18 | 函数平移导致次微分定义域相应平移 |
| 4.3 | 18–19 | \(\max\{g(x),-\sqrt t\}\) 的凸性及次微分定义域公式 |
| 5.1 | 19–20 | 缩放集合 SDF 给出的加法上界 |
| 5.2 | 20–21 | 一般底集生成锥的法锥与极锥公式 |
| 5.3 | 21 | 一般底集生成锥的边界投影候选分段 |

### 引理（Lemma）

| 编号 | PDF 页 | 内容 |
|---|---:|---|
| 3.1 | 7 | 开球切点与切半空间外一点的凸包，会把切点变为内部点；用于法向几何证明 |
| 3.2 | 7–8 | 两个不同内部点若共用一个最近边界点，则三点共线 |
| 3.3 | 9 | 骨架为空当且仅当是半空间或没有内部 |
| 3.4 | 10 | 内部点由唯一 V(C) 起点沿某最近边界方向分解 |
| 3.5 | 11 | 单位法向半空间 \(D_H=\langle a,x\rangle-b\) |
| 3.6 | 12 | 闭凸集为全部支撑半空间之交 |
| 3.7 | 15 | \(\operatorname{conv}\operatorname{dom}f=\operatorname{dom}\operatorname{conv}f\) |
| 3.8 | 15 | proper 凸函数的定义域在取函数闭包前后相对内部相同 |

### 推论（Corollary）

| 编号 | PDF 页 | 内容 |
|---|---:|---|
| 3.1 | 6 | 适当正则条件下，次微分和的值域与次微分值域之和具有相同 SDF |
| 3.2 | 6 | proper 凸集 SDF 的 Eikonal 性质几乎处处成立 |
| 3.3 | 10 | 闭凸集的边界投影参数单调 |
| 3.4 | 12 | SDF 共轭为球面限制支撑函数的闭凸化 |
| 3.5 | 13 | nearly convex 集合的 SDF 与其闭包、相对内部相同，且凸 |
| 3.6 | 13 | 归一化多面体约束在内部给出有限 max 公式 |
| 3.7 | 15 | 有内部时，SDF 极小点必属于骨架 |
| 3.8 | 16 | 紧凸集 SDF 次微分值域为整个单位球 |
| 3.9 | 16 | 闭凸锥 SDF 次微分值域为极锥单位方向的凸包 |
| 3.10 | 16 | 无内部时 SDF 次微分值域等于普通距离值域公式 |
| 3.11 | 17 | 外部所有梯度恰为实际法锥值域的单位向量 |
| 4.1 | 18 | 可逆线性变换下次微分定义域按同一变换变化 |
| 4.2 | 18 | 次微分定义域可以进行任意可逆仿射变换 |
| 4.3 | 18 | 次微分定义域类对等距变换封闭 |
| 4.4 | 19 | 任意朝向的半空间删除其边界一块闭凸片的相对内部后，仍属于次微分定义域类 |

## 12. 原文笔误及使用时的修正

以下基于 PDF 页面图像核对及定义直接验证，属于本解读的勘误说明；没有据此断言其他结论错误。

1. **PDF 第 17 页，期刊 p.767，例 3.7：半空间边界投影的正负号。** 原文写为 \(x+a(\langle a,x\rangle-b)\)，正确应为

   \[
   Q_H(x)=\{x-a(\langle a,x\rangle-b)\},\quad\|a\|=1.
   \]

   代回即可检查 \(\langle a,Q_H(x)\rangle=b\)。例如 H 为 \(x_1\le0\)，(1,0) 应投影到 (0,0)，原文加号却给出 (2,0)。

2. **PDF 第 13 页，期刊 p.763，推论 3.6 的证明：多面体边界写法。** 原文写 \(\partial C=\bigcap_i\{x:a_i^Tx=b_i\}\)，一般不对。正确为

   \[
   \partial C=\bigcup_i\big(C\cap\{x:a_i^Tx=b_i\}\big).
   \]

   正方形四个面所属超平面的交甚至为空。该处错误不改变内部 SDF 的 max 公式，可从“到补集各开放半空间的距离取最小值”正确证明。

3. **PDF 第 21 页，期刊 p.771，例 5.2：二阶锥 Q 的分段遗漏。** 原文对 \(r>0\) 的投影只列 \(-r<t\le r\)，漏掉 \(t>r\) 的非轴内部点。实际相同公式 \(\tfrac12(r+t)(z/r,1)\) 对所有 \(r>0,t>-r\) 成立。次微分分段也应补上负轴非原点，避免除零；第 8.4 节已给完整表达。

4. **PDF 第 5 页，期刊 p.755，命题 3.2 证明中的复合函数下标。** 一般相似变换的正确法则为 \(D_{\phi(S)}=\alpha D_S\circ\phi^{-1}\)。原文证明在右侧重复写成 \(D_{\phi(S)}\)，命题列出的具体变换结论正确。

另外，“几乎处处单位梯度”不是“处处光滑”；“SDF 凸”不是“SDF 处处可微”；“凸函数”不是“次微分定义域必凸”；边界投影 Q 与集合投影 P 在内部不同。这四处概念区分比记忆公式更重要。

## 13. 应用与理论结论的对应

下列计算场景是对论文结论的解释或直接推论，并非论文分别提出的新算法。

| 场景 | 需要的结论 | 具体计算 |
|---|---|---|
| 点与障碍物的间距、内部深度 | 定义、法向分解、次微分 | 符号判断内外；绝对值为边界距离；可微时还原接触点 |
| 带半径的小球与闭障碍物 C | SDF 定义 | 额外假设外部球心 p：严格不相交需 \(D_C(p)>r\)，接触对应等号；若要求球位于 C 内则 \(D_C(p)\le-r\) |
| 可行域的安全裕量 | 凸性、多面体内部式 | \(D_C\le-\varepsilon\) 编码距边界至少 \(\varepsilon\) 的内部点；凸 C 时仍为凸约束 |
| 最大内接球 | 定理 3.9 | 最小化 D，负最小值为半径；最优点满足最近边界点凸包条件 |
| 水平集、形状分析 | \(\{D=0\}=\partial C\)、Eikonal 性质 | 用零水平集表达边界，用单位梯度描述距离校准 |
| 二分类几何分数 | 内负外正与距离大小 | 在已有区域模型下，符号提供标签，大小表示几何边界余量；不等同概率置信度 |
| 凸约束投影 | 式 (11)、正交锥、二阶锥公式 | 多面体用二次规划；标准锥可直接计算 |
| 凸分析与单调算子反例 | 第 4 章构造 | 有效定义域凸，但有限次梯度存在域可非凸 |

若只需要掌握本文理论主线，优先读：定义与 Fact 1 → 定理 3.1 → 定理 3.3 → 定理 3.7 → 定理 3.8 → 定理 3.9 → 定理 4.1 → 定理 5.1。其余结果解释投影与次微分的精细结构。

<!-- sdf-terminology-answer:start -->
## 附录：第 2 节 Terminology and Notation 详细解释

> 对应问题：详细解释一下 PDF 的 2 Terminology and Notation 中的内容，例如提到的 convex 的定义等等。

第 2 节的作用，是为后面讨论 SDF 建立一套共同语言。它把**集合的几何、函数的凸性，以及多值算子的性质**放到同一个框架里。

下面对应原文 PDF 第 3—4 页，即期刊页码 753—754；定义按原文解释，例子和推导是补充说明。原文为 Honglin Luo、Xianfu Wang、Brett Lukens 的 *Variational Analysis on the Signed Distance Functions*，Journal of Optimization Theory and Applications 180 (2019), 751–774，DOI：10.1007/s10957-018-1414-2。

### 1．基本空间与符号：实向量空间、内积、范数、单位球

\(\mathbb R^n\) 是所有 \(n\) 维实向量组成的空间。例如：

\[
x=(x_1,\ldots,x_n)\in\mathbb R^n.
\]

二维对应平面，三维对应通常的三维空间。

论文使用标准欧氏内积：

\[
\langle x,y\rangle=\sum_{i=1}^n x_i y_i,
\]

以及由它定义的欧氏范数：

\[
\|x\|=\sqrt{\langle x,x\rangle}
=\sqrt{x_1^2+\cdots+x_n^2}.
\]

例如 \(x=(3,4)\)，则 \(\|x\|=5\)。两点间的距离就是 \(\|x-y\|\)。

内积还反映方向：

\[
\langle x,y\rangle=\|x\|\|y\|\cos\theta.
\]

因此，内积为正表示夹角小于 \(90^\circ\)，为零表示垂直，为负表示夹角大于 \(90^\circ\)。后面的“法锥”和“单调性”都依靠这个几何意义。

论文还使用以下记号：

| 记号 | 定义 | 二维中的形状 |
|---|---|---|
| \(\mathbb B\) | \(\{x:\|x\|\le1\}\) | 包含边界的单位圆盘 |
| \(\mathbb B^\circ\) | \(\{x:\|x\|<1\}\) | 不包含边界的单位圆盘 |
| \(\mathbb S\) | \(\{x:\|x\|=1\}\) | 单位圆周 |
| \(\mathbb R_+\) | \([0,+\infty)\) | 非负实数 |
| \(\mathbb R_{++}\) | \((0,+\infty)\) | 正实数 |

注意区分：\(\mathbb S\) 是固定的单位球面，而普通字母 \(S\) 可以表示任意集合。

### 2．线段与凸集：convex 的定义究竟是什么意思

给定两个点 \(x,y\)，论文定义连接它们的闭线段：

\[
[x,y]
=\{(1-\lambda)x+\lambda y:0\le\lambda\le1\}.
\]

其中：

- \(\lambda=0\) 时得到 \(x\)；
- \(\lambda=1\) 时得到 \(y\)；
- \(\lambda=\frac12\) 时得到中点；
- 其他值对应线段上的其他位置。

例如 \(x=(0,0),y=(2,4)\)，则

\[
(1-\lambda)x+\lambda y=(2\lambda,4\lambda).
\]

当 \(\lambda\) 从 0 变到 1 时，该点沿线段从 \(x\) 移动到 \(y\)。

论文用

\[
]x,y[=[x,y]\setminus\{x,y\}
\]

表示去掉两个端点的开线段。

**凸集的定义是：**

\[
\boxed{
x,y\in C,\quad0\le\lambda\le1
\quad\Longrightarrow\quad
(1-\lambda)x+\lambda y\in C.
}
\]

直观地说：**在集合内任取两个点，连接它们的整条线段都留在集合内。**

这里有两个“任意”：任意两个点、任意比例。只要找到一组点和比例不满足条件，就能证明集合不凸。

例如：

- 区间 \([-1,1]\) 是凸集。
- 圆盘、矩形、半空间是凸集。
- 圆周不是凸集：\((1,0)\) 与 \((-1,0)\) 在圆周上，但中点 \((0,0)\) 不在。
- 圆环不是凸集：连接两侧的点，线段会穿过中间的空洞。

以单位圆盘为例，若 \(\|x\|\le1,\|y\|\le1\)，根据三角不等式：

\[
\|(1-\lambda)x+\lambda y\|
\le(1-\lambda)\|x\|+\lambda\|y\|
\le1.
\]

因此整条线段仍在圆盘内，这就从定义证明了圆盘的凸性。

**凸性与是否包含边界、是否有界，是不同性质。**

| 集合 | 凸？ | 闭？ | 有界？ |
|---|---:|---:|---:|
| 开圆盘 | 是 | 否 | 是 |
| 闭圆盘 | 是 | 是 | 是 |
| 半空间 | 是 | 是 | 否 |
| 圆周 | 否 | 是 | 是 |

原文凸集定义中的 \(\subset\) 应理解为“包含于”，不要求是真包含。

论文还称集合 \(C\) 为 **proper**，如果

\[
\varnothing\ne C\ne\mathbb R^n.
\]

它只表示集合既非空，也不是整个空间；不代表集合凸、闭或有界。空集和整个空间本身也满足凸集定义，只是被 proper 条件排除了。

### 3．补集、闭包、内部与边界

这些概念用来准确表达“内部、外部、边界”。

**补集：**

\[
S^c=\mathbb R^n\setminus S.
\]

即所有不属于 \(S\) 的点。

**闭包：**

\[
\operatorname{cl}S
\]

是 \(S\) 加上它所有极限点得到的集合。可以理解为“把缺失的边界补齐”。

例如：

\[
S=(0,1)
\quad\Longrightarrow\quad
\operatorname{cl}S=[0,1].
\]

虽然 \(0\notin S\)，但 \(1/k\in S\)（取 \(k\ge2\)）且 \(1/k\to0\)，所以 \(0\) 属于闭包。

**内部：**

点 \(x\) 是内部点，是指存在一个以 \(x\) 为中心的小开球，完全包含在 \(S\) 中：

\[
x+\varepsilon\mathbb B^\circ\subseteq S
\quad\text{对某个 }\varepsilon>0.
\]

所有内部点组成 \(\operatorname{int}S\)。

**边界：**

论文定义

\[
\boxed{
\operatorname{bdry}S
=\operatorname{cl}S\cap\operatorname{cl}(S^c).
}
\]

边界点的特点是：无论在它周围取多小的邻域，都能遇到集合中的点，也能遇到补集中的点。

对 \(S=(0,1)\)：

\[
\operatorname{int}S=(0,1),\qquad
\operatorname{cl}S=[0,1],\qquad
\operatorname{bdry}S=\{0,1\}.
\]

“内部”取决于所在空间。例如二维中的一条直线没有二维内部，因为任何小圆盘都会离开直线。这正是后文某些集合的 SDF 不会出现负值的原因。

### 4．扩展实数、有效定义域和 proper 函数

论文允许函数取无穷值：

\[
\overline{\mathbb R}
=\mathbb R\cup\{-\infty,+\infty\}.
\]

为什么要允许 \(+\infty\)？因为这样可以把约束写进函数中。

例如：

\[
f(x)=
\begin{cases}
x^2,&x\ge0,\\
+\infty,&x<0.
\end{cases}
\]

最小化这个函数，就相当于求解

\[
\min_{x\ge0}x^2.
\]

负数位置的 \(+\infty\) 表示这些输入不允许被采用。

论文定义**有效定义域**：

\[
\operatorname{dom}f
=\{x:f(x)<+\infty\}.
\]

对上面的例子：

\[
\operatorname{dom}f=[0,+\infty).
\]

虽然函数形式上写成 \(f:\mathbb R\to\overline{\mathbb R}\)，但真正允许有限代价的输入只有这个集合。

函数 \(f\) 是 **proper**，需要：

\[
f(x)>-\infty\quad\forall x,
\qquad
\operatorname{dom}f\ne\varnothing.
\]

即：

1. 任何位置都不取 \(-\infty\)；
2. 至少有一个位置取有限值。

所以 \(f(x)=x^2\) 是 proper；恒等于 \(+\infty\) 的函数不是 proper。

**集合 proper 与函数 proper 是两套不同定义。** 比如 \(x^2\) 是 proper 函数，尽管其有效定义域就是整个 \(\mathbb R\)。

论文还约定：

\[
\sup\varnothing=-\infty,\qquad
\inf\varnothing=+\infty.
\]

这样距离和支撑函数即使遇到空集，也有统一的表达式。例如到空集的距离为 \(+\infty\)。

### 5．上图集 epigraph：如何把函数问题变成集合问题

论文定义：

\[
\boxed{
\operatorname{epi}f
=\{(x,r)\in\mathbb R^n\times\mathbb R:r\ge f(x)\}.
}
\]

它由函数图像上的点，以及所有位于图像上方的点组成。

例如 \(f(x)=x^2\)：

\[
\operatorname{epi}f=\{(x,r):r\ge x^2\}.
\]

于是：

- \((2,4)\) 在上图集边界上；
- \((2,5)\) 在上图集内部；
- \((2,3)\) 不在上图集中。

需要区分：

\[
\text{函数图像：}\quad r=f(x),
\]

\[
\text{上图集：}\quad r\ge f(x).
\]

引入上图集后，“函数凸不凸”“函数闭不闭”都能转化为一个集合的性质。

### 6．凸函数：为什么论文用“上图集是凸集”来定义

原文定义：

\[
\boxed{
f\text{ 是凸函数}
\iff
\operatorname{epi}f\text{ 是凸集}.
}
\]

它等价于常见的不等式：

\[
\boxed{
f((1-\lambda)x+\lambda y)
\le
(1-\lambda)f(x)+\lambda f(y),
\quad0\le\lambda\le1.
}
\]

这里可以先对 \(x,y\in\operatorname{dom}f\) 理解该不等式；在这个范围内，proper 函数的端点值均为有限实数。

几何含义是：**连接函数图像上任意两点的线段，都位于函数图像上方或图像上。**

为什么两个定义等价？

取上图集中的两点：

\[
(x,f(x)),\qquad(y,f(y)).
\]

它们的凸组合为

\[
\bigl((1-\lambda)x+\lambda y,\,
(1-\lambda)f(x)+\lambda f(y)\bigr).
\]

如果上图集凸，这个点仍在上图集中，因此它的高度必须不低于对应函数值，恰好就得到上述凸性不等式。

例如 \(f(x)=x^2\)：

\[
(1-\lambda)x^2+\lambda y^2
-\bigl((1-\lambda)x+\lambda y\bigr)^2
=
\lambda(1-\lambda)(x-y)^2\ge0.
\]

所以 \(x^2\) 凸。

相反，\(f(x)=-x^2\) 不凸。取 \(x=-1,y=1,\lambda=\tfrac12\)，凸性要求

\[
f(0)\le\frac{f(-1)+f(1)}2,
\]

但实际是 \(0\le-1\)，不成立。

几个容易混淆的地方：

- 凸函数不必可微，例如 \(|x|\)。
- 凸函数不必单调，例如 \(x^2\)。
- 仿射函数也凸，例如 \(2x+1\)。
- 凸函数的**上图集**凸，其函数图像一般不凸。
- proper 凸函数的有效定义域一定是凸集。

### 7．闭函数、下半连续、函数闭包与函数凸包

**闭函数 closed**

原文定义：

\[
f\text{ 闭}
\iff
\operatorname{epi}f\text{ 闭}.
\]

这不是说函数图像闭，也不是只检查 \(\operatorname{dom}f\) 是否闭。

它等价于**下半连续**，缩写为 lsc：

\[
\boxed{
x_k\to x
\quad\Longrightarrow\quad
f(x)\le\liminf_{k\to\infty}f(x_k).
}
\]

直观上：函数在某点的值，不能突然高于邻近点能够逼近的函数值。

例如：

\[
f(0)=1,\qquad f(x)=0\quad(x\ne0).
\]

取 \(x_k=1/k\)，则

\[
f(0)=1>0=\liminf f(x_k),
\]

所以它不下半连续。

反过来，若

\[
g(0)=0,\qquad g(x)=1\quad(x\ne0),
\]

则它下半连续。下半连续允许这种向下的突降。

**函数的闭包 \(\operatorname{cl}f\)**

论文通过上图集定义：

\[
\operatorname{epi}(\operatorname{cl}f)
=\operatorname{cl}(\operatorname{epi}f).
\]

即先补齐上图集缺失的极限点，再把结果视为一个函数的上图集。

例如令

\[
f(x)=
\begin{cases}
0,&0<x<1,\\
+\infty,&\text{其他}.
\end{cases}
\]

取闭包后变成：

\[
(\operatorname{cl}f)(x)=
\begin{cases}
0,&0\le x\le1,\\
+\infty,&\text{其他}.
\end{cases}
\]

端点处缺失的有限值被补上了。

**函数的凸包 \(\operatorname{conv}f\)**

先回顾集合的凸包：\(\operatorname{conv}A\) 是包含 \(A\) 的最小凸集，也就是所有有限凸组合构成的集合。

论文定义：

\[
(\operatorname{conv}f)(x)
=
\inf\{r:(x,r)\in\operatorname{conv}(\operatorname{epi}f)\}.
\]

操作是：**先对上图集取凸包，再取其下边界。**

例如：

\[
f(-1)=1,\qquad f(0)=2,\qquad f(1)=1,
\]

其他位置取 \(+\infty\)。则

\[
(\operatorname{conv}f)(x)=
\begin{cases}
1,&-1\le x\le1,\\
+\infty,&\text{其他}.
\end{cases}
\]

因为连接 \((-1,1)\) 和 \((1,1)\) 后，会把中间高度为 1 的线段加入上图集。

所以，函数的凸化通常是**从下方降低函数，并可能填充定义域的空隙**。

### 8．Fenchel 次微分：用支撑斜率代替导数

原文定义：

\[
\boxed{
\partial f(x)
=
\left\{
v:
f(z)\ge f(x)+\langle v,z-x\rangle
\quad\forall z
\right\}.
}
\]

其中 \(x\in\operatorname{dom}f\)；如果 \(x\notin\operatorname{dom}f\)，则规定

\[
\partial f(x)=\varnothing.
\]

固定 \(x,v\) 后，

\[
f(x)+\langle v,z-x\rangle
\]

是一个经过 \((x,f(x))\) 的仿射函数。

因此，\(v\) 是次梯度的意思是：**它给出一张接触函数、且在整个空间内都不超过函数的支撑平面。**

对于可微凸函数：

\[
\partial f(x)=\{\nabla f(x)\}.
\]

但次微分也能处理不可微点。

以 \(f(x)=|x|\) 为例，在 \(x=0\)：

\[
v\in\partial f(0)
\iff
|z|\ge vz\quad\forall z.
\]

- 对 \(z>0\)，得到 \(v\le1\)；
- 对 \(z<0\)，得到 \(v\ge-1\)。

所以：

\[
\boxed{\partial |x|(0)=[-1,1].}
\]

完整表达为：

\[
\partial |x|=
\begin{cases}
\{-1\},&x<0,\\
[-1,1],&x=0,\\
\{1\},&x>0.
\end{cases}
\]

尖点处没有唯一导数，却存在一整组支撑斜率。

还要区分：

\[
\operatorname{dom}f
\quad\text{和}\quad
\operatorname{dom}\partial f.
\]

前者表示函数有有限值，后者表示至少存在一个有限次梯度。

例如：

\[
f(x)=
\begin{cases}
-\sqrt{x},&x\ge0,\\
+\infty,&x<0.
\end{cases}
\]

它是 proper 闭凸函数，且 \(0\in\operatorname{dom}f\)。但若 \(v\in\partial f(0)\)，必须有

\[
-\sqrt z\ge vz\quad\forall z>0,
\]

即

\[
v\le-\frac1{\sqrt z}.
\]

当 \(z\to0^+\) 时右侧趋于 \(-\infty\)，没有有限 \(v\) 能满足，所以

\[
\partial f(0)=\varnothing.
\]

**函数有值，不保证该点有次梯度。** 这正是文章第 4 章研究次微分定义域的基础。

### 9．多值映射、定义域、值域、图像与逆

原文写：

\[
A:\mathbb R^n\rightrightarrows\mathbb R^n.
\]

双箭头表示 \(A(x)\) 是一个集合，可以为空，也可以含一个或多个向量。次微分就是多值映射。

它的定义域、值域、图像分别是：

\[
\operatorname{dom}A=\{x:A(x)\ne\varnothing\},
\]

\[
\operatorname{ran}A=\bigcup_x A(x),
\]

\[
\operatorname{gr}A=\{(x,u):u\in A(x)\}.
\]

以 \(A=\partial|\cdot|\) 为例：

\[
\operatorname{dom}A=\mathbb R,
\qquad
\operatorname{ran}A=[-1,1].
\]

图像包含三部分：负半轴对应高度 \(-1\)，正半轴对应高度 \(1\)，以及 \(x=0\) 上从 \(-1\) 到 \(1\) 的竖线段。

多值逆定义为：

\[
\boxed{A^{-1}(u)=\{x:u\in A(x)\}.}
\]

它问的是：“哪些输入能够产生输出 \(u\)？”

例如：

\[
A^{-1}(0)=\{0\},\qquad
A^{-1}(1)=[0,+\infty),\qquad
A^{-1}(2)=\varnothing.
\]

这里的逆不是数值倒数，也不要求输出唯一。

### 10．monotone、maximally monotone 与 paramonotone

这三个定义都针对算子的图像。

**单调 monotone**

\[
\boxed{
\langle x-y,u-v\rangle\ge0,
\quad
u\in A(x),\ v\in A(y).
}
\]

在一维中就是：

\[
(x-y)(u-v)\ge0.
\]

也就是说，输入增大时，允许的输出不能发生相反方向的变化。这推广了普通函数的单调不减。

在高维中，它约束的是两个变化向量的内积，而不是要求每个分量分别递增。

例如 \(A(x)=x\)：

\[
\langle x-y,A(x)-A(y)\rangle=\|x-y\|^2\ge0,
\]

所以单调。\(A(x)=-x\) 则不单调。

凸函数的次微分为何单调？若

\[
u\in\partial f(x),\qquad v\in\partial f(y),
\]

则

\[
f(y)\ge f(x)+\langle u,y-x\rangle,
\]

\[
f(x)\ge f(y)+\langle v,x-y\rangle.
\]

两式相加，就得到

\[
\langle x-y,u-v\rangle\ge0.
\]

**极大单调 maximally monotone**

它表示：算子已经单调，并且**无法再向图像中添加新的输入输出对而保持单调**。

“极大”针对图像的包含关系，不是函数值大，也不是斜率大。

例如定义：

\[
A(x)=
\begin{cases}
\{-1\},&x<0,\\
\{0\},&x=0,\\
\{1\},&x>0.
\end{cases}
\]

它单调，但不是极大单调。因为可以把 \(A(0)\) 扩充成 \([-1,1]\)，而不破坏单调性。

扩充后恰好得到 \(\partial|x|\)。这个完整算子是极大单调的。

**参数单调 paramonotone**

它先要求单调，再要求：

\[
\langle x-y,u-v\rangle=0
\quad\Longrightarrow\quad
u\in A(y),\ v\in A(x).
\]

即：**如果单调性不等式恰好取等号，两点的输出可以交换。**

它不涉及额外的“参数”，名称不能按字面理解。

一个非平凡例子是：

\[
f(x_1,x_2)=|x_1|.
\]

对 \(x=(0,0)\)、\(y=(0,1)\)，两处次微分都是

\[
[-1,1]\times\{0\}.
\]

选择 \(u=(1,0),v=(-1,0)\)，则

\[
\langle x-y,u-v\rangle
=\langle(0,-1),(2,0)\rangle=0.
\]

两个输出确实都能在另一点取得。

但单调性本身不保证参数单调。例如旋转算子

\[
A(x_1,x_2)=(-x_2,x_1)
\]

满足

\[
\langle h,A(h)\rangle=0,
\]

所以单调。可是取 \(x=(1,0),y=(0,0)\)，两点输出分别为 \((0,1)\) 和 \((0,0)\)，不能交换，因此不参数单调。

极大单调关注“图像是否还能扩充”；参数单调关注“等号情况下能否交换输出”。两者研究不同问题。

### 11．与集合相关的三个函数：示性函数、支撑函数、距离函数

这三个函数从不同角度表达同一个集合。

**示性函数 indicator function**

\[
\boxed{
\iota_S(x)=
\begin{cases}
0,&x\in S,\\
+\infty,&x\notin S.
\end{cases}}
\]

这里是 \(0/+\infty\)，不是概率论中常见的 \(0/1\)。

它用于把约束并入目标函数：

\[
\min_{x\in S}g(x)
\quad\Longleftrightarrow\quad
\min_x\bigl(g(x)+\iota_S(x)\bigr).
\]

\(S\) 凸，当且仅当 \(\iota_S\) 凸。

**支撑函数 support function**

\[
\boxed{
\sigma_S(v)=\sup_{s\in S}\langle v,s\rangle.
}
\]

固定方向 \(v\)，它衡量集合在这个方向上的最远投影位置。

例如 \(S=[-1,1]\)：

\[
\sigma_S(v)=|v|.
\]

因为 \(v>0\) 时最远点为 \(1\)，\(v<0\) 时最远点为 \(-1\)。

对单位球：

\[
\sigma_{\mathbb B}(v)=\|v\|.
\]

这来自 Cauchy–Schwarz 不等式：

\[
\langle v,s\rangle\le\|v\|\|s\|\le\|v\|,
\]

且可以在 \(v\ne0\) 时选择与 \(v\) 同方向的单位向量达到等号；\(v=0\) 时等式显然成立。

支撑函数是方向性的几何描述，不是点到集合的距离。

**距离函数 distance function**

\[
\boxed{d_S(x)=\inf_{s\in S}\|x-s\|.}
\]

例如 \(S=[-1,1]\)：

\[
d_S(x)=\max\{|x|-1,0\}.
\]

距离为零只保证 \(x\in\operatorname{cl}S\)，未必属于 \(S\)。例如 \(S=(0,1)\)，则 \(d_S(0)=0\)，但 \(0\notin S\)。

同样的区间 \(S=[-1,1]\)，三种表达分别是：

\[
\iota_S(x)=0\text{ 或 }+\infty,
\qquad
\sigma_S(v)=|v|,
\qquad
d_S(x)=\max\{|x|-1,0\}.
\]

它们分别表达：**是否允许、方向范围、最近距离。**

### 12．法锥 normal cone：为什么是“锥”，为什么朝外

按照论文定义：

\[
N_S(x)=
\begin{cases}
\{v:\langle v,s-x\rangle\le0,\ \forall s\in S\},
&x\in\operatorname{cl}S,\\
\varnothing,&x\notin\operatorname{cl}S.
\end{cases}
\]

重点看：

\[
\langle v,s-x\rangle\le0.
\]

\(s-x\) 是从 \(x\) 指向集合中其他点的向量。上述条件说明：当相关向量非零时，\(v\) 与这些向量的夹角都不小于 \(90^\circ\)。

因此，非零 \(v\) 指向集合的外侧，并给出一张经过 \(x\) 的支撑超平面。

为什么叫“锥”？因为如果 \(v\) 满足条件，那么任意 \(\alpha\ge0\) 都满足：

\[
\langle\alpha v,s-x\rangle
=\alpha\langle v,s-x\rangle\le0.
\]

所以法锥收集的是一整组方向及其所有非负倍数，而不是只有单位法向量。

例如 \(C=[-1,1]\)：

\[
N_C(x)=
\begin{cases}
\{0\},&-1<x<1,\\
[0,+\infty),&x=1,\\
(-\infty,0],&x=-1,\\
\varnothing,&x\notin C.
\end{cases}
\]

在右端点 \(1\)，所有向右的向量都是外法向；在左端点则相反。

对单位圆盘，若 \(\|q\|=1\)：

\[
N_C(q)=\{\alpha q:\alpha\ge0\}.
\]

对正方形右上角 \((1,1)\)，法锥是整个第一象限：

\[
N_C(1,1)=\{(a,b):a\ge0,b\ge0\}.
\]

这解释了为什么光滑边界通常只有一个单位外法向，而角点可以有很多。

对于闭凸集，还有一个重要联系：

\[
\boxed{N_C=\partial\iota_C.}
\]

法锥可以看作示性函数的次微分。因此，集合的几何约束与函数的支撑斜率，本质上可以用同一套工具处理。

### 13．普通投影、边界投影与 SDF

普通投影定义为：

\[
P_S(x)=\{s\in S:\|x-s\|=d_S(x)\}.
\]

距离 \(d_S(x)\) 只是一个数，投影 \(P_S(x)\) 则记录**哪些点达到这个距离**。

这里使用集合，是因为最近点可能不唯一，也可能根本不存在。例如 \(S=(0,1)\)、\(x=0\)，距离为零，但没有集合中的点达到零距离，所以 \(P_S(0)=\varnothing\)。

在有限维欧氏空间中，非空闭集的投影一定存在；闭凸集的普通投影还具有唯一性。

论文随后定义边界投影：

\[
\boxed{Q_S(x)=P_{\operatorname{bdry}S}(x).}
\]

它寻找最近的边界点。对于 proper 集合，边界闭且非空，所以这个投影始终存在。

最后，SDF 定义为：

\[
\boxed{
D_S(x)=
\begin{cases}
d_S(x),&x\notin S,\\
-d_{S^c}(x),&x\in S.
\end{cases}}
\]

普通距离在集合内部为零，SDF 则以负数表达内部深度。

仍以 \(C=[-1,1]\) 为例：

\[
D_C(x)=|x|-1.
\]

| 查询点 | 普通距离 \(d_C(x)\) | SDF \(D_C(x)\) | 普通投影 \(P_C(x)\) | 边界投影 \(Q_C(x)\) |
|---|---:|---:|---|---|
| \(x=2\) | \(1\) | \(1\) | \(\{1\}\) | \(\{1\}\) |
| \(x=0.3\) | \(0\) | \(-0.7\) | \(\{0.3\}\) | \(\{1\}\) |
| \(x=0\) | \(0\) | \(-1\) | \(\{0\}\) | \(\{-1,1\}\) |
| \(x=1\) | \(0\) | \(0\) | \(\{1\}\) | \(\{1\}\) |

中心点 \(0\) 有两个最近边界点，这与 \(|x|-1\) 在零点不可微相对应；但它的次微分仍存在：

\[
\partial D_C(0)=[-1,1].
\]

这个例子把第 2 节的主要概念连接起来：**集合的边界投影可能多值，而 SDF 的次微分正好能够表达这种几何上的不唯一性。**
<!-- sdf-terminology-answer:end -->

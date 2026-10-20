# versor-java

面向计算机图形学的 Java 几何代数（Clifford 代数）库。

API 语义参考 `versor`（C++，`/home/wangwenjie/Documents/代码库/versor`：
`Multivector` 的 `gp / op / ip / dual / reverse`，`Gen::rot / trs / mot / log`），
采用 Java 的运行期签名 + 位掩码实现，不照搬其编译期模板类型系统。
数值用 PyPI 的 `galgebra` 与 Wolfram `CliffordAlgebra` 交叉验证。

## 范围

- 任意 `Cl(p,q,r)`，生成元 ≤ 14 个（$2^n$ 系数）；`Signature.cl(p,q)` 直接给出
  前 p 个 +1、后 q 个 −1 的度量；n≤10 用预计算符号表，n=11–14 按需计算符号
  （每对更慢，避免 1 GiB 级符号表）
- 多向量：几何积、外积、左/右缩并、回归积、交换子/反交换子、
  reverse/conjugate/involution、grade 投影、dual/undual、标量积、奇 versor 反射、
  多向量除法、范数、逆（versor 快速路径与带均衡化/残差检查的任意可逆元素通用求逆）、
  指数（缩放平方和）、`isFinite/isScalar/isHomogeneous/isBlade/isVersor` 校验
- `MultivectorBatch`：按 blade 连续存储的批量运算（gp / gpRight / op / ip /
  sandwich / dual / undual / reverse / 算术），避免逐元素对象分配；
  `sandwich(versor, out)` 复用输出缓冲，`slice(from,to)` 支持分块并行
- `SparseMultivector`：COO（mask + value）表示，低密度积用哈希累加，
  高密度退回稠密累加器；适合 k≪2^n 的高维输入
- `MultivectorMatrix`：多向量矩阵乘（i-k-j，复用稠密乘积累加器）、转置与
  Gauss-Jordan 求逆
- `E3`：向量、转子、旋转、slerp、log
- `P3`（PGA3）：点/面/线、translator、motor、变换、平面/直线求交、
  线方向/矩/锚点、点到线/面距离、最近点、平面投影、motor 插值
- `Conformal`（Cl(4,1)）：点/球/平面/圆、球坐标解码、点到球面/平面距离与最近点
- `Motion`：twist 定义的刚体运动，解析的位置/速度/加速度/jerk
- `Field` / `Calculus`：多重向量场，∇F、方向导数、散度、旋度
- `Manifold`：参数化曲线/曲面，切刃与法向
- 不包含：versor 的绘图/App 层，CGA 的 conic/cyclide/tangent 等高级形态

## 架构

| 层 | 类型 | 表示与职责 |
|---|---|---|
| 代数 | `Signature` | 度量、blade 位掩码、grade、乘积符号；所有操作数必须共享同一实例 |
| 稠密 | `Multivector` | 按 mask 索引的 `double[2^n]`；`gp` 自适应选择稠密循环或非零支撑对 |
| 稀疏 | `SparseMultivector` | COO `(mask, value)`；按规模在哈希/稠密累加器间切换 |
| 批量 | `MultivectorBatch` | SoA 按 blade 连续存储；`sandwich(versor,out)` 复用输出，`slice` 分块 |
| 矩阵 | `MultivectorMatrix` | i-k-j 多向量矩阵乘、转置、Gauss-Jordan 求逆 |
| 模型 | `E3` / `P3` / `Conformal` | 向量、转子、motor、点/面/线及 fast path |
| 运动 | `Motion` / `Field` / `Calculus` / `Manifold` | twist 运动学、场微分、参数化流形 |

不变量：blade 按整数 mask 排序；零度量生成元的平方为 0；n≤10 使用惰性符号表，
n=11–14 按需计算符号。批量和矩阵操作要求形状与签名一致。

内核分派与统一分析见 [`docs/KERNEL_SCHEDULING.md`](docs/KERNEL_SCHEDULING.md)。

## 约定（务必对齐）

PGA3 基序为 `e0, e1, e2, e3`，其中 `e0² = 0`，`e1² = e2² = e3² = 1`。

| 元素 | 表示 |
|---|---|
| 点 | `P(x,y,z) = e123 - z e012 + y e013 - x e023`（三矢量） |
| 面 `n·x + d = 0` | `π = nx e1 + ny e2 + nz e3 + d e0`（矢量） |
| 平移 | `T = 1 - (tx e01 + ty e02 + tz e03)/2`，`T P ~T` 平移 `+t` |
| 旋转 | `R = cos(θ/2) - sin(θ/2)(ax e23 - ay e13 + az e12)` |
| motor | `M = T R`（先旋转，后平移），`M P ~M` |
| 线 | 两点 `p1.vee(p2)`；两面 `π1.op(π2)`（对偶表示，二矢量） |
| 面-线交点 | `P3.meetPlaneLine(π, L) = π.op(L)` |

## 用法

```java
import vsr.*;

// E3: 绕 z 轴 90°，把 e1 转到 e2
Multivector r = E3.rotor(0, 0, 1, Math.PI / 2);
Multivector img = E3.rotate(E3.vec(1, 0, 0), r);   // -> e2

// 转子插值
Multivector mid = Gen.slerp(E3.rotor(0,0,1,0), r, 0.5);

// PGA3: 刚体变换
Multivector m = P3.motor(
        P3.rotor(0, 0, 1, Math.PI / 2),
        P3.translator(1, 2, 3));
double[] q = P3.pointCoords(P3.apply(m, P3.point(1, 0, 0)));

// 刚体运动插值（旋转 slerp + 平移线性）
Multivector half = P3.slerp(Multivector.scalar(P3.SIG, 1), m, 0.5);

// 求交：平面 z=2 与 z 轴
Multivector line = P3.meetPlanes(
        P3.plane(1, 0, 0, 0), P3.plane(0, 1, 0, 0));
double[] hit = P3.pointCoords(
        P3.meetPlaneLine(P3.plane(0, 0, 1, -2), line)); // (0,0,2)

// 几何查询：点到平行 z 轴的直线的距离与最近点
Multivector zAxis = P3.line(P3.point(1, 2, 0), P3.point(1, 2, 1));
double[] closest = P3.pointCoords(
        P3.closestPointOnLine(zAxis, P3.point(0, 0, 0))); // (1,2,0)
double distance = P3.distancePointLine(zAxis, P3.point(0, 0, 0));

// CGA：点到球面的距离与最近点
Multivector cgaSphere = Conformal.sphere(1, 1, 1, 2);
double[] sc = Conformal.sphereCoords(cgaSphere);   // {1, 1, 1, 2}
double surfaceDistance = Conformal.distancePointSphere(cgaSphere, Conformal.point(1, 1, 1));
double[] onSphere = Conformal.pointCoords(
        Conformal.closestPointOnSphere(cgaSphere, Conformal.point(1, 1, 1))); // (3,1,1)
```

示例：`examples/CubeDemo.java` 把一个立方体的 8 个顶点用 motor 变换并打印。

### C++ versor 语义对齐

```java
Multivector a = E3.vec(1, 2, 3);
Multivector b = E3.vec(4, 5, 6);

double dot = a.scalarProduct(b);          // 32
Multivector q = a.divide(b);              // a * b^-1（b 为 versor/blade）
Multivector back = q.mul(b);              // 回到 a

Multivector d = E3.vec(1, 0, 0).dual();
Multivector same = d.undual();            // 补集对偶的逆

Multivector reflected = E3.vec(0, 1, 0).reflect(E3.vec(1, 0, 0)); // e2 不变
```

`undual()` 是 Java 补集 `dual()` 的逆，不等同于 C++ 的 Hodge `A * I`；
退化度量下不要直接套用 C++ 的 `dual/undual` 符号。

### 可读写法

Java 没有运算符重载，除 C++ 对应的 `gp/op/ip/vee/sandwich` 外，提供一组读法更直接的别名：

| 可读写法 | 等价写法 | 语义 |
|---|---|---|
| `a.mul(b)` | `a.gp(b)` | 几何积 |
| `a.mul(s)` | `a.scale(s)` | 标量乘法 |
| `a.wedge(b)` | `a.op(b)` | 外积 |
| `a.contract(b)` | `a.ip(b)` | 左缩并 |
| `a.regressive(b)` | `a.vee(b)` | 回归积 |
| `versor.apply(x)` | `versor.sandwich(x)` | 偶 versor 变换 |

```java
Multivector rotor = E3.rotor(0, 0, 1, Math.PI / 2);
Multivector image = rotor.apply(E3.vec(1, 0, 0));   // e2
Multivector area = E3.vec(1, 0, 0).wedge(E3.vec(0, 1, 0)); // e12
```

## 运动学：速度、加速度、jerk

`Motion` 用 twist（二矢量生成元 `W`）定义 `M(t) = exp(tW)`。对点
`P(t) = M p ~M`，导数全部解析：

```text
P'   = W P - P W
P''  = W^2 P + P W^2 - 2 W P W
P''' = W^3 P - 3 W^2 P W + 3 W P W^2 - P W^3
```

```java
Motion m = Motion.screw(0, 0, 1, 2.0, 0, 0, 1); // 绕 z 角速度 2 + 沿 z 速度 1
double[] p = m.position(0.5, 1, 0, 0);
double[] v = m.velocity(0.5, 1, 0, 0);
double[] a = m.acceleration(0.5, 1, 0, 0);
double[] j = m.jerk(0.5, 1, 0, 0);

Motion spin = Motion.rotation(0, 0, 1, Math.PI);   // 纯旋转
Motion drift = Motion.translation(1, 0, 0);        // 纯平移
```

示例：`examples/KinematicsDemo.java`。

## Blaze3D / Minecraft 集成（可选）

`src/blaze3d/java` 是 JOML 适配层（Blaze3D 的 `PoseStack` 消费
`Quaternionf` / `Matrix4f` / `Vector3f`）。核心保持零依赖，适配层单独编译：

```bash
./build-blaze3d.sh     # 自动查找 PrismLauncher / ~/.m2 里的 joml，或 JOML_JAR=/path/joml.jar
```

```java
// PoseStack：先旋转后平移，与 M = T R 的顺序一致
poseStack.mulPose(Blaze3D.quaternion(P3.rotationOf(motor)));
Vector3f t = Blaze3D.translation(motor);
poseStack.translate(t.x, t.y, t.z);

// 或直接取矩阵 / 变换顶点
Matrix4f m = Blaze3D.matrix(motor);
Vector3f p = Blaze3D.transform(motor, new Vector3f(1, 0, 0));
Multivector back = Blaze3D.motor(m);              // 矩阵 -> P3 motor
Quaternionf q = Blaze3D.quaternion(rotor);        // E3/P3 rotor -> 四元数
```

`Matrix4f` 按列主序（JOML 的 `m<col><row>`）；绕 z 轴 90° 时 `m01 = 1`、`m10 = -1`。

## 场与流形（几何微积分）

`Calculus.vectorDerivative` 直接计算

\[ \nabla F = e_1\partial_1 F + e_2\partial_2 F + e_3\partial_3 F \]

标量部分是散度，二矢量部分取对偶即旋度：

```java
Field swirl = (x, y, z) -> E3.vec(-y, x, 0);
double div = Calculus.divergence(swirl, 0, 0, 0);   // 0
Multivector curl = Calculus.curl(swirl, 0, 0, 0);   // 2 e3
Multivector dv = Calculus.directional(scalarField, 1, 2, 3, 1, 0, 0);

Manifold sphere = Manifold.surface(u -> new double[]{
        2 * Math.cos(u[0]) * Math.cos(u[1]),
        2 * Math.sin(u[0]) * Math.cos(u[1]),
        2 * Math.sin(u[1])});
Multivector n = sphere.normal(0, 0);                // 单位法向 e1
```

`Manifold.tangentBlade` 对曲线返回切向量、对曲面返回切二矢量；
`normal` 只在曲面上定义，取切二矢量的对偶并归一化。
`jacobian` 给出曲线速度 / 曲面面积元，`Calculus.flux` 用中点法算
\(\int F\cdot n\,dS\)，这是粒子按面积加权采样的基础。

## 参考资料

| 来源 | 用途 | 采用情况 |
|---|---|---|
| Alan Macdonald, *Linear and Geometric Algebra* (2010)，中文译本 | 转子 \(R=e^{-B\theta/2}\)、旋转/反射、复合与插值（第 7 章） | 用于 `E3`/`P3`/`Gen` 的符号与手性校验 |
| Alan Macdonald, *Vector and Geometric Calculus* (2012) | 微分、方向导数、向量导数、线/面积分（第 3、5、7、9 章） | `Calculus` 的 \(\nabla F\)、散度、旋度；扫掠与面积加权采样 |
| [jeremyong/klein](https://github.com/jeremyong/klein)（MIT，已归档） | PGA 标量参考、测试黄金值、SIMD lane 布局 | 参考符号表/展开思路；**未复制代码**；SIMD 仅作可选原型，不进入模组 |
| [MozillaFiredoge/Aerodynamics4MC-Core](https://github.com/MozillaFiredoge/Aerodynamics4MC-Core)（MIT） | 薄入口、服务端权威、量化同步、粒子风场 | 采用架构模式；**未复制代码** |
| Prodanov, [*Algorithmic computation of multivector inverses...*](https://arxiv.org/abs/2308.02291)（2023） | Faddeev–LeVerrier–Souriau 多向量逆/特征多项式 | 采用“先判定可逆、再求逆”的思路；实现仍用数值线性代数 |
| Higham, [*Iterative Refinement Implies Numerical Stability for Gaussian Elimination*](https://doi.org/10.1090/S0025-5718-1980-0572859-4)（1980） | 迭代精化与后向误差 | `inverseGeneral` 的残差校正与检查准则 |
| [LAPACK Users' Guide, *How to Measure Errors*](https://netlib.org/lapack/lug/node75.html) | 条件数、后向误差、`RCOND` | 条件数估计与误差度量 |
| Garamon, [*A Geometric Algebra Library Generator*](https://doi.org/10.1007/s00006-019-0987-7) | 代码生成/展开的性能上限 | 固定公式展开路线，未复制代码 |
| [*Java Vector API: Benchmarking and Performance Analysis*](https://doi.org/10.1145/3578360.3580265)（2023） | Vector API 的布局与收益边界 | SIMD 原型仅用于批量 SoA 内核 |

本地书籍路径：`/home/wangwenjie/Documents/论文资料/books/`（不进仓库）。

## 性能

`bench.sh` 是零依赖基准（维度扫描、批量扩展、分配量、矩阵/高维压力、C1/C2 对比，
检测到 JOML 时附带四元数/矩阵横向对比）。`bench-jmh.sh` 下载 JMH 1.37 并跑
带 `-prof gc` 的复核；`render` 模式覆盖端到端、并发与机器人求解：

```bash
./bench.sh              # 全量
./bench.sh sweep        # 维度扫描
./bench.sh sparse       # 高维稀疏 n=8–14
./bench.sh support      # 稠密数组上的稀疏支撑快速路径
./bench.sh render       # 端到端/并发/机器人（需要 JOML）
./bench.sh robust       # 结构化/混合密度/策略交叉/符号表
./bench.sh batchscale   # 分配 vs 复用、线程与块大小
./bench.sh multicore    # gp/op/batch/scalar 多核扩展
./bench.sh fuzz         # 稠密-稀疏/条件数/退化/精度（需要 JOML）
./bench.sh simd         # 可选 jdk.incubator.vector 批量 sandwich 原型
./bench.sh gabm         # 输出 ga-benchmark 可读 JSON（需要 JOML）
./run-gabm.sh           # 生成 + 校验；设置 GA_BENCHMARK_DIR 可直接跑 gabm.py
# 生成的 JSON 含 mean/stddev 聚合，可直接交给 ga-benchmark 的 gabm.py
./bench.sh jmh 'JmhDense' -prof gc
./bench.sh jmh 'JmhFast' -t 1
./bench-all.sh          # 全套矩阵（测试、sweep、sparse、support、simd、robust、render、批量、fuzz、多核、JMH、ga-benchmark）
BENCH_JAVA_OPTS="-XX:TieredStopAtLevel=1" ./bench.sh jit
```

完整结果、环境与方法见 [`PERFORMANCE.md`](PERFORMANCE.md)；排版版本见
[`PERFORMANCE.pdf`](PERFORMANCE.pdf)（R 生成表格与图，本地 XeLaTeX 编译）：

```bash
cd report
python3 import_log.py /tmp/full-new.log
Rscript report2.R                 # 需要 R + xelatex + 思源宋体/黑体
```

## 批量与矩阵

```java
// 批量:元素 i 只依赖输入 i,系数按 blade 连续存放
MultivectorBatch batch = MultivectorBatch.of(
        P3.point(1, 0, 0), P3.point(0, 1, 0), P3.point(0, 0, 1));
MultivectorBatch moved = batch.sandwich(motor);   // motor * p * ~motor
Multivector p1 = moved.get(1);

// 多向量矩阵:乘法按几何积,求逆用左侧行变换 + 通用主元求逆
MultivectorMatrix m = new MultivectorMatrix(P3.SIG, 3, 3);
m.set(0, 0, P3.motor(P3.rotor(0, 0, 1, 0.5), P3.translator(1, 0, 0)));
MultivectorMatrix inv = m.inverse();              // m.inverse() 要求方阵且可逆

// 通用求逆:均衡化 + 两次迭代精化 + 后向误差检查
Multivector a = Multivector.scalar(E3.SIG, 1.3)
        .add(E3.vec(0.4, -0.2, 0.7)).add(E3.bivector(0.1, 0.3, -0.2));
Multivector invA = a.inverseGeneral();             // 相对残差 > 1e-12 时抛异常
Multivector.GeneralInverse info = a.inverseGeneralResult(); // 不抛，读 residual/condition
```

`MultivectorBatch` 的运算复杂度是 `O(4^n · count)`;`MultivectorMatrix.inverse()`
对每个主元调用 `inverseGeneral()`(系数空间 `O(8^n)`),适合中小规模压力测试。

`inverseGeneral()` 对系数矩阵做行列均衡化，按最大行/列尺度选主元，再用
`A x = 1` 的残差做至多两次迭代精化。默认要求相对后向误差 ≤ `1e-12` 且条件数 ≤ `1e12`；
需要保留病态结果时用 `inverseGeneralResult()` 读取 `residual`、`condition`，
或传自定义容差 `inverseGeneral(tol)`。条件数很大时即使残差合格，前向误差仍会
被放大；机器人 IK/物理仿真应优先使用阻尼最小二乘或 SVD 伪逆。

## 构建与测试

```bash
./build.sh        # javac + 自带断言测试（无需 Maven/Gradle）
```

`build.sh` 会编译 `src/main/java`，然后运行 `vsr.GaTest`（454 项断言）。
`build-blaze3d.sh` 额外编译 JOML 适配层并运行 `vsr.blaze3d.Blaze3DTest`（11 项断言）。
示例：

```bash
javac -cp out -d out examples/CubeDemo.java && java -cp out CubeDemo
```

## 验证

- galgebra（PyPI）与 Wolfram `CliffordAlgebra` 的黄金值：E3 多向量乘积、
  转子系数、`R e1 ~R` / `R e2 ~R`
- 随机结合律、`(AB)~ = B~ A~`、`A∧A` 的反对称性
- PGA：点坐标往返、平移/旋转/复合、逆、面-面/面-线求交、点连线
- 通用逆：尺度不变、残差/条件数可见、近奇异结果可拒绝
- 转子：`R ~R = 1`、`exp(log R) = R`、slerp 中点
- motor：`T R` 作用顺序、逆、`translationOf`/`rotationOf`、插值中点
- 运动学：旋转/平移/螺旋的解析位置、速度、加速度、jerk，并与位置差分、
  加速度差分交叉验证
- 几何微积分：旋量场散度/旋度/∇F、标量场梯度、方向导数；平面/球面法向、
  螺旋线切向；面积元与通量积分（球面通量 = 4πR³）

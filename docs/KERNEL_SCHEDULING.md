# 内核统一与调度分析

本文分析 `versor-java` 自身的运算内核、当前分派点，以及本轮为减少重复和中间
分配所做的统一。目标不是引入外部线性代数库，而是让已有的稠密、稀疏、批量和
矩阵路径共享同一套系数级内核，并在正确的粒度上选择路径。

## 1. 当前层次

| 层 | 实现 | 存储 | 主要分派 |
|---|---|---|---|
| 符号 | `Signature` | 度量/位掩码/惰性符号表 | n≤10 表，n=11–14 按需计算 |
| 稠密 | `Multivector` | `double[2^n]` | `gp` 按非零支撑选稀疏或全扫描；`op/ip/rc` 枚举兼容子掩码 |
| 稀疏 | `SparseMultivector` | COO 排序 mask/value | 低密度 open-addressing，高密度稠密累加器 |
| 批量 | `MultivectorBatch` | SoA `blades[mask][i]` | E3/PGA3 展开 sandwich；一般路径两次批量积 |
| 矩阵 | `MultivectorMatrix` | `Multivector[rows][cols]` | 多向量几何积 + Gauss-Jordan |

关键不变量：同一操作中的所有操作数共享一个 `Signature`；blade 由整数 mask
排序；`Signature` 负责乘积符号和 dual 符号，避免稠密、批量和稀疏各自维护
一套符号逻辑。

## 2. 现有的调度点

### 2.1 稠密 `gp`

`Multivector.gp` 先统计两个操作数的非零系数数 $n_a,n_b$。当

$$
4 n_a n_b < 2^n \cdot 2^n
$$

时只遍历非零支撑，否则走全扫描。这个阈值把稀疏散列/索引构建成本与
$4^n$ 次乘法做了粗略比较。实际选择已经由 `PERFORMANCE.md` 的稀疏基准验证。

### 2.2 外积、缩并与右缩并

`op`、`ip`、`rightContraction` 不再枚举全部 $4^n$ 对，而是枚举与左操作数
兼容的子掩码：

- $A\wedge B$：只枚举与 $A$ 不相交的 blade；
- $A\rfloor B$：只枚举包含 $A$ 的 blade；
- $A\lfloor B$：只枚举 $A$ 的子掩码。

这把 n=10 的外积从 $4^{10}$ 对降到 $3^{10}$ 对。

### 2.3 稀疏 `gp`

`SparseMultivector.gp` 根据 $n_a n_b$ 选择哈希累加器或稠密累加器。阈值
`hashFactor` 在 n≤10 时为 4，在更高维为 16，反映符号表查询成本的变化。

### 2.4 批量 sandwich

`MultivectorBatch.sandwich` 对 E3 向量和 PGA3 点使用与标量路径同形的展开
公式；一般 versor 使用两次批量积：

$$
M X \tilde M = (M X)\tilde M.
$$

这样比三重循环更省算术，代价是中间 batch。输出复用版本
`sandwich(versor, out)` 对一般路径保留三重展开，避免分配中间 batch。

### 2.5 矩阵

`MultivectorMatrix.multiply` 以前按 `r-c-k` 顺序对每个结果元素累加
`Multivector.add(gp(...))`，每个 $k$ 产生一个临时多向量。本轮改为：

1. `i-k-j` 顺序，左因子每个 `(row,k)` 只读取一次；
2. 每个结果元素只分配一个 `double[2^n]` 累加器；
3. 复用 `Multivector.gpAccumulate`，不再为每个 $k$ 构造 `Multivector`。

## 3. 本轮的统一

### 3.1 共享稠密乘积累加器

`Multivector.gp` 的循环体移到包内方法

```java
void gpAccumulate(Multivector o, double[] out)
```

`gp` 只负责分配结果数组，矩阵乘法负责提供结果累加器。这样稠密乘积只有一份
实现，同时把“输出分配”与“系数累加”解耦。`gp` 的稀疏支撑选择仍在累加器
内部，因此矩阵路径也自动获得同一分派。

### 3.2 矩阵循环调度

`MultivectorMatrix.multiply` 的新结构：

```text
for r:
  for k:
    left = A[r][k]
    if left == 0: continue
    for c:
      left.gpAccumulate(B[k][c], acc[r][c])
```

相比原来的 `r-c-k`，它把同一个左因子保持在寄存器/缓存中，并让 `B[k]` 行
连续访问。结果矩阵每个 cell 只构造一次 `Multivector`。

矩阵消元也从 `sub(factor.gp(cell))` 改为“复制目标 cell + 以 `-1` 为 scale 的
`gpAccumulate`”，去掉一次中间多向量分配。主元归一化仍保留已确认的
`inverseGeneral` 路径。

### 3.3 批量支撑扫描

`MultivectorBatch.gp` 和 `product` 在进入主要乘积循环前预先计算 batch 中哪些
blade 非零；一般 sandwich 的输出复用路径也只扫描一次。这样内层循环不再对
每个固定操作数重复 `hasNonzero`。

`product` 的固定左/右操作数分派统一成同一张 `active[]` 支撑表；E3/PGA3
展开公式仍保留独立实现，因为它们的主要收益来自消除通用乘积的索引和符号
计算，而不是共享控制流。

## 4. 本轮测量

单线程微基准，Linux/JDK 25，非隔离环境：

| 路径 | 旧 | 新 | 变化 |
|---|---:|---:|---:|
| P3 3×3 多向量矩阵乘（20k 次） | 4820.2 ns | 3744.7 ns | **1.29×** |
| Cl(4) 批量一般 sandwich，输出复用（10k 点，300 次） | 1.04 ms（分配+复制） | 0.93 ms | **1.12×** |
| Cl(4) 批量一般 sandwich，分配路径 | 0.80 ms（两次批量积） | 1.00 ms（三重展开） | **0.80×** |

最后一行说明：**不能为了统一而强行让所有入口走三重展开**。分配路径保留
两次批量积；只有调用者显式提供 `out` 时才使用输出复用三重展开。这个分支
差异是性能必要的，不是重复。

## 5. 仍需统一/优化的地方

1. **E3/PGA3 展开公式**：`Multivector` 与 `MultivectorBatch` 各有一份。
   可以生成一份符号项表（term plan），标量路径和批量路径都按项表执行；
   但直接共享控制流会引入虚调用/边界检查，必须先测。当前保留两份展开，
   用测试锁定一致性。
2. **矩阵求逆**：Gauss-Jordan 的归一化和消元仍对每个 cell 调用 `gp`/`sub`。
   可以用 `gpAccumulate` + 负因子累加消除中间 `Multivector`，再做 3×3
   专用展开。
3. **稀疏 sandwich**：`PERFORMANCE.md` 已指出高密度中间结果会退化；下一步
   应直接展开三重积并按中间支撑剪枝，而不是继续调 hash/dense 阈值。
4. **阈值集中**：稠密 `4 n_a n_b`、稀疏 `hashFactor` 和批量 `N` 阈值目前
   分散在各自类中。可以收进一个包内 `KernelPolicy`，但只有出现第三个调用方
   时才值得抽象；当前先保留局部常量并加注释。
5. **并行调度**：批量已经通过 `slice` 支持外部并行；库内不默认创建线程池，
   避免和 Minecraft 服务端/客户端线程竞争。并行阈值应按 chunk 内元素数和
   内存带宽决定，而不是简单按线程数。

## 6. 验证

```bash
JAVA_HOME=~/.jdks/jbrsdk-25.0.3 ./gradlew gaTest
JAVA_HOME=~/.jdks/jbrsdk-25.0.3 ./gradlew combatTest
JAVA_HOME=~/.jdks/jbrsdk-25.0.3 ./gradlew weTest
```

`gaTest` 覆盖稠密/稀疏/批量/矩阵路径；矩阵乘法和批量输出复用的一致性由
现有测试和新增的微基准交叉检查。

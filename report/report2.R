#!/usr/bin/env Rscript
# versor-java 2.0 report generator. Run after bench-all.sh imports its CSVs.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(script_dir)
dir.create("figures2", showWarnings = FALSE)

sweep <- read.csv("data/sweep.csv")
batch <- read.csv("data/batch.csv")
stress <- read.csv("data/stress.csv")
jit <- read.csv("data/jit.csv")
joml <- read.csv("data/joml.csv")
bigsparse <- read.csv("data/big-sparse.csv")
robust <- read.csv("data/robust.csv")
threshold <- read.csv("data/threshold.csv")
sigtab <- read.csv("data/table.csv")
render <- read.csv("data/render.csv")
robot <- read.csv("data/robot.csv")
batchscale <- read.csv("data/batchscale.csv")
fuzz <- read.csv("data/fuzz.csv")
multicore <- read.csv("data/multicore.csv")
jmh <- read.csv("data/jmh.csv")
submask <- read.csv("data/submask.csv")
support <- if (file.exists("data/support.csv")) read.csv("data/support.csv") else NULL
gabm <- if (file.exists("data/gabm.csv")) read.csv("data/gabm.csv") else NULL
meta <- read.csv("data/meta2.csv")
meta_val <- function(k) meta$value[meta$key == k]

fmt <- function(x) ifelse(x >= 1000, formatC(round(x), format = "d", big.mark = ","), sprintf("%.1f", x))
tex <- function(x) gsub("_", "\\\\_", x)
fixnl <- function(x) {
  x <- gsub("\\n", "\n", x, fixed = TRUE)
  v <- strsplit(x, "\n", fixed = TRUE)[[1]]
  v <- sub("\\\\+$", "\\\\\\\\", v)
  paste(v, collapse = "\n")
}
CJK <- "Source Han Sans CN"

gp10 <- sweep$ns[sweep$n == 10 & sweep$op == "gp"]
pair10 <- gp10 / 4^10
op10 <- sweep$ns[sweep$n == 10 & sweep$op == "op"]
sub10 <- submask[submask$n == 10 & submask$op == "op", ]
sp10 <- bigsparse$ns[bigsparse$n == 10 & bigsparse$k == 10 & bigsparse$op == "gp"]
sp14 <- bigsparse[bigsparse$n == 14 & bigsparse$k == 512 & bigsparse$op == "gp", ]
sand14 <- bigsparse[bigsparse$n == 14 & bigsparse$k == 512 & bigsparse$op == "sandwich", ]
rga <- render[render$variant == "ga", ]
rb <- render[render$variant == "batch", ]
rmx <- render[render$variant == "matrix", ]
rad <- render[render$variant == "adapter", ]
alloc1m <- batchscale[batchscale$section == "reuse" & batchscale$a == "1000000" & batchscale$b == "alloc", ]
reuse1m <- batchscale[batchscale$section == "reuse" & batchscale$a == "1000000" & batchscale$b == "reuse", ]
mc20 <- function(w) multicore[multicore$workload == w & multicore$threads == 20, ]

summary <- paste0(
"\\begin{itemize}\n",
"\\item 稠密几何积仍是唯一 $4^n$ 热点：n=10 为 ", sprintf("%.2f", gp10 / 1e6), " ms，每对约 ",
sprintf("%.2f", pair10), " ns；输出驱动与 FMA 变体实测更慢（n=10 为 1.47--1.50 ms），保留现有循环。\n",
"\\item 核心 \\texttt{gp} 增加自适应稀疏支撑：PGA motor 乘点 ", fmt(support$ns[support$kind == "P3"]),
" ns；Cl(10) 仅 10 个非零系数时 ", sprintf("%.1f", support$ns[support$kind == "sparse10" & support$k == 10] / 1000),
" $\\mu$s，而稠密 ", sprintf("%.2f", gp10 / 1e6), " ms。\n",
"\\item 外积/缩并改为子掩码枚举（$4^n\\to3^n$）：n=10 外积 ", sprintf("%.1f", op10 / 1000),
" $\\mu$s，原 ", sprintf("%.0f", sub10$before / 1000), " $\\mu$s（", sprintf("%.1f", sub10$before / sub10$after), "$\\times$）。\n",
"\\item 稀疏 COO：n=10、k=10 的 \\texttt{gp} ", sprintf("%.2f", sp10 / 1000), " $\\mu$s，约 ",
sprintf("%.0f", gp10 / sp10), "$\\times$ 于稠密；n=14、k=512 的 \\texttt{gp}/\\texttt{sandwich} 为 ",
sprintf("%.0f", sp14$ns / 1e6), "/", sprintf("%.0f", sand14$ns / 1e6), " ms。\n",
"\\item 批量输出复用把 1M 点从 ", sprintf("%.1f", alloc1m$ns), " 降到 ",
sprintf("%.1f", reuse1m$ns), " ns/点（分配 128 $\\to$ 0 B/点）；4 线程达到带宽平台。\n",
"\\item 端到端 1M 点：GA ", sprintf("%.1f", rga$ns_per_point), "、批量 ",
sprintf("%.1f", rb$ns_per_point), "、缓存 Matrix4f ", sprintf("%.2f", rmx$ns_per_point),
"、adapter ", sprintf("%.0f", rad$ns_per_point), " ns/点。\n",
"\\item 多核 20 线程：gp ", sprintf("%.2f", mc20("gp_n10")$speedup), "$\\times$、op ",
sprintf("%.2f", mc20("op_n10")$speedup), "$\\times$、批量 ", sprintf("%.2f", mc20("batch_sandwich")$speedup),
"$\\times$、标量 ", sprintf("%.2f", mc20("scalar_sandwich")$speedup), "$\\times$。\n",
"\\item 断言：", meta_val("checks"), " 项通过；稠密/稀疏逐位一致，子掩码与暴力参考一致；",
"固定主元阈值在近奇异时可能静默返回不精确逆。\n",
"\\end{itemize}\n")

# ---------------------------------------------------------------- figures
cairo_pdf("figures2/fig-gp.pdf", width = 6.2, height = 4.0, family = CJK)
gp <- sweep[sweep$op == "gp", ]; gp <- gp[order(gp$n), ]
plot(gp$n, gp$ns, log = "xy", type = "b", pch = 19, lwd = 2, col = "#1f5fa8",
     xlab = "代数维度 n", ylab = "单次 gp (ns, 对数轴)", main = "稠密 gp 与 4^n 参考")
lines(gp$n, gp$ns[1] * 4^(gp$n - gp$n[1]), lty = 2, col = "grey40")
legend("topleft", c("实测", "4^n"), lty = c(1, 2), pch = c(19, NA), col = c("#1f5fa8", "grey40"), bty = "n")
dev.off()

cairo_pdf("figures2/fig-submask.pdf", width = 6.2, height = 3.6, family = CJK)
ops <- c("op", "ip", "rightContraction"); cols <- c("#1f5fa8", "#b45309", "#3f7f5f")
plot(NA, xlim = range(submask$n), ylim = c(1, max(submask$before / submask$after) * 1.15),
     xlab = "代数维度 n", ylab = "加速比", main = "子掩码枚举：4^n → 3^n")
for (i in seq_along(ops)) { d <- submask[submask$op == ops[i], ]; d <- d[order(d$n), ]
  lines(d$n, d$before / d$after, type = "b", pch = 19, col = cols[i], lwd = 2) }
legend("topleft", ops, pch = 19, col = cols, bty = "n")
dev.off()

cairo_pdf("figures2/fig-threshold.pdf", width = 6.2, height = 3.8, family = CJK)
th <- threshold[threshold$n == 14, ]
plot(th$k, th$hash_ns, log = "xy", type = "b", pch = 19, col = "#b45309", lwd = 2,
     xlab = "非零数 k (对数轴)", ylab = "gp (ns, 对数轴)", main = "n=14：哈希 vs 稠密累加")
lines(th$k, th$dense_ns, type = "b", pch = 17, col = "#1f5fa8", lwd = 2)
lines(th$k, th$auto_ns, lty = 2, col = "grey30", lwd = 1.5)
legend("topleft", c("哈希", "稠密", "自动"), pch = c(19, 17, NA), lty = c(1, 1, 2),
       col = c("#b45309", "#1f5fa8", "grey30"), bty = "n")
dev.off()

cairo_pdf("figures2/fig-batch.pdf", width = 6.2, height = 3.8, family = CJK)
lp <- batch$loop_ns / batch$n; bp <- batch$batch_ns / batch$n
plot(batch$n, lp, log = "xy", type = "b", pch = 19, lwd = 2, col = "#1f5fa8",
     ylim = range(c(lp, bp)), xlab = "批量大小 N", ylab = "每点 ns (对数轴)",
     main = "逐元素 vs SoA 批量")
lines(batch$n, bp, type = "b", pch = 17, lwd = 2, col = "#b45309")
legend("topright", c("Multivector[]", "MultivectorBatch"), pch = c(19, 17),
       col = c("#1f5fa8", "#b45309"), bg = "white", bty = "o", box.col = NA, cex = 0.9)
dev.off()

cairo_pdf("figures2/fig-reuse.pdf", width = 6.2, height = 3.6, family = CJK)
r1 <- batchscale[batchscale$section == "reuse" & batchscale$a == "1000000", ]
par(mar = c(4, 7, 3, 1))
barplot(rev(r1$ns), names.arg = rev(r1$b), horiz = TRUE, las = 1, col = "#3f7f5f",
        xlab = "ns/点", main = "100 万点：分配 vs 输出复用")
dev.off()

cairo_pdf("figures2/fig-multicore.pdf", width = 6.2, height = 3.8, family = CJK)
workloads <- c("gp_n10", "op_n10", "batch_sandwich", "scalar_sandwich")
cols <- c("#1f5fa8", "#b45309", "#3f7f5f", "#7c3aed")
plot(NA, xlim = range(multicore$threads), ylim = c(0, 20), xlab = "线程数", ylab = "加速比",
     main = "多核扩展（20 逻辑核）")
abline(a = 0, b = 1, lty = 2, col = "grey50")
for (i in seq_along(workloads)) { d <- multicore[multicore$workload == workloads[i], ]
  d <- d[order(d$threads), ]
  lines(d$threads, d$speedup, type = "b", pch = 19, col = cols[i], lwd = 2) }
legend("topleft", workloads, pch = 19, col = cols, bty = "n", cex = 0.9)
dev.off()

cairo_pdf("figures2/fig-joml.pdf", width = 6.2, height = 3.8, family = CJK)
sel <- c("E3,rotate", "JOML,rotateAlloc", "P3,applyPoint", "JOML,matrixTransformAlloc")
labels <- c("E3 rotate", "JOML quat", "P3 motor×point", "JOML matrix")
vals <- sapply(sel, function(s) { p <- strsplit(s, ",", fixed = TRUE)[[1]];
  joml$ns[joml$path == p[1] & joml$op == p[2]][1] })
par(mar = c(4, 10, 3, 1))
barplot(vals, names.arg = labels, horiz = TRUE, las = 1, log = "x", col = "#3f7f5f",
        xlab = "ns/op（对数轴）", main = "GA 与 JOML")
dev.off()

cairo_pdf("figures2/fig-cond.pdf", width = 6.2, height = 3.6, family = CJK)
cd <- fuzz[fuzz$cat == "cond", ]
plot(as.numeric(cd$b), pmax(cd$value, 1e-30), log = "xy", type = "b", pch = 19,
     col = "#b45309", lwd = 2, xlab = "扰动 eps（对数轴）", ylab = "残差（对数轴）",
     main = "inverseGeneral 近奇异残差")
dev.off()

# ----------------------------------------------------------------- tables
sweep_table <- function() {
  ops <- c("gp", "op", "ip", "rightContraction", "commutator", "reverse", "dual",
           "norm2", "unit", "inverse", "sandwich", "inverseGeneral", "exp")
  dims <- sort(unique(sweep$n))
  head <- paste0("\\textbf{运算} & ", paste0("\\textbf{n=", dims, "}", collapse = " & "), " \\\\")
  rows <- vapply(ops, function(op) {
    cells <- vapply(dims, function(n) { v <- sweep$ns[sweep$op == op & sweep$n == n]
      if (!length(v)) "--" else fmt(v[1]) }, "")
    paste0("\\texttt{", tex(op), "} & ", paste(cells, collapse = " & "), " \\\\") }, "")
  paste0("\\resizebox{\\textwidth}{!}{%\\n{\\scriptsize\\setlength{\\tabcolsep}{3pt}\\n",
         "\\begin{tabular}{l", strrep("r", length(dims)), "}\\n\\toprule\\n", head,
         "\\n\\midrule\\n", paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}}}\\n")
}

fast_table <- function() {
  gp3 <- sweep$ns[sweep$n == 3 & sweep$op == "gp"]
  gp4 <- sweep$ns[sweep$n == 4 & sweep$op == "gp"]
  e3 <- joml$ns[joml$path == "E3" & joml$op == "rotate"]
  p3 <- joml$ns[joml$path == "P3" & joml$op == "applyPoint"]
  rows <- c(sprintf("E3 rotate & %s & 通用 n=3 gp: %s \\\\", fmt(e3), fmt(gp3)),
            sprintf("P3 motor×point & %s & 通用 n=4 gp: %s \\\\", fmt(p3), fmt(gp4)))
  paste0("\\begin{tabular}{lll}\\n\\toprule\\n路径 & ns/op & 对照 \\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

support_table <- function() {
  if (is.null(support)) return("")
  sp <- support[support$kind == "sparse10", ]
  rows <- paste0(sp$k, " & ", fmt(sp$ns), " \\\\")
  paste0("\\begin{tabular}{rr}\\n\\toprule\\nCl(10) 非零数 k & gp ns \\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

submask_table <- function() {
  d <- submask[submask$n %in% c(3, 6, 8, 10), ]
  rows <- paste0(d$n, " & \\texttt{", tex(d$op), "} & ", fmt(d$before), " & ",
                 fmt(d$after), " & ", sprintf("%.2f$\\times$", d$before / d$after), " \\\\")
  paste0("\\begin{tabular}{rlrrr}\\n\\toprule\\nn & 运算 & 原 ns & 新 ns & 提升 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

bigsparse_table <- function() {
  pick <- bigsparse[bigsparse$op == "gp" & bigsparse$n %in% c(10, 12, 14)
                    & bigsparse$k %in% c(8, 32, 512), ]
  rows <- vapply(seq_len(nrow(pick)), function(i) { r <- pick[i, ]
    s <- bigsparse$ns[bigsparse$n == r$n & bigsparse$k == r$k & bigsparse$op == "sandwich"]
    sprintf("%d & %d & %s & %s \\\\", r$n, r$k, fmt(r$ns),
            if (length(s)) fmt(s[1]) else "--") }, "")
  paste0("\\begin{tabular}{rrrr}\\n\\toprule\\nn & k & gp ns & sandwich ns \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

threshold_table <- function() {
  d <- threshold[threshold$n %in% c(10, 14) & threshold$k %in% c(2, 8, 16, 32, 128, 512), ]
  rows <- paste0(d$n, " & ", d$k, " & ", fmt(d$hash_ns), " & ", fmt(d$dense_ns),
                 " & ", fmt(d$auto_ns), " & ", d$choice, " & ", d$correct, " \\\\")
  paste0("\\begin{tabular}{rrrrrll}\\n\\toprule\\nn & k & hash & dense & auto & 选择 & 判定 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

sigtab_table <- function() {
  d <- sigtab[sigtab$threads == 1 & sigtab$n <= 10, ]
  rows <- paste0(d$n, " & ", sprintf("%.1f", d$first_us), " & ", sprintf("%.1f", d$rebuild_us),
                 " & ", ifelse(d$bytes == 0, "--", sprintf("%.1f MiB", d$bytes / 1048576)), " \\\\")
  paste0("\\begin{tabular}{rrrl}\\n\\toprule\\nn & 首次乘积 us & 重建 us & 表内存 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

robust_table <- function() {
  d <- robust[robust$k %in% c(8, 32)
              & robust$kind %in% c("uniform", "versor", "pgaPoint")
              & robust$n %in% c(10, 14), ]
  keys <- unique(d[, c("n", "kind", "k")])
  rows <- vapply(seq_len(nrow(keys)), function(i) {
    r <- keys[i, ]
    gp <- d$ns[d$n == r$n & d$kind == r$kind & d$k == r$k & d$op == "gp"]
    sw <- d$ns[d$n == r$n & d$kind == r$kind & d$k == r$k & d$op == "sandwich"]
    sprintf("%d & \\texttt{%s} & %d & %s & %s \\\\", r$n, tex(r$kind), r$k,
            if (length(gp)) fmt(gp[1]) else "--",
            if (length(sw)) fmt(sw[1]) else "--")
  }, "")
  paste0("\\begin{tabular}{rlrrr}\\n\\toprule\\nn & 结构 & k & gp ns & sandwich ns \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

stress_table <- function() {
  d <- stress[order(stress$n), ]
  f <- function(v) ifelse(v <= 0, "--", fmt(v))
  rows <- paste0(d$n, " & ", f(d$dense_sandwich_ns), " & ",
                 f(d$inverse_general_ns), " & ", f(d$matrix3_inverse_ns), " \\\\")
  paste0("\\begin{tabular}{rrrr}\n\\toprule\nn & 稠密 sandwich & inverseGeneral & 3x3 矩阵逆 \\\\\\n\\midrule\n",
         paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

jmh_table <- function() {
  d <- jmh[jmh$bench == "JmhDense" & jmh$threads == 1
            & jmh$op %in% c("gp", "sandwich", "inverse")
            & jmh$n %in% c(3, 4, 6, 8, 10), ]
  keys <- sort(unique(d$n))
  rows <- vapply(keys, function(n) {
    g <- d$ns[d$n == n & d$op == "gp"]; sw <- d$ns[d$n == n & d$op == "sandwich"]
    iv <- d$ns[d$n == n & d$op == "inverse"]
    sprintf("%d & %s & %s & %s \\\\", n,
            if (length(g)) fmt(g[1]) else "--",
            if (length(sw)) fmt(sw[1]) else "--",
            if (length(iv)) fmt(iv[1]) else "--")
  }, "")
  f <- jmh[jmh$bench == "JmhFast" & jmh$n == 0, ]
  fast <- c("e3Rotate", "p3Transform", "loopTransform", "batchTransform", "inverseGeneral")
  frows <- vapply(fast, function(op) {
    a <- f$ns[f$op == op & f$threads == 1]; b <- f$ns[f$op == op & f$threads == 4]
    sprintf("\\texttt{%s} & %s & %s \\\\", tex(op),
            if (length(a)) fmt(a[1]) else "--", if (length(b)) fmt(b[1]) else "--")
  }, "")
  paste0("\\begin{minipage}{0.56\\textwidth}\\centering\\small\n",
         "\\begin{tabular}{rrrr}\n\\toprule\nn & gp & sandwich & inverse \\\\\\n\\midrule\n",
         paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n\\end{minipage}\n",
         "\\begin{minipage}{0.40\\textwidth}\\centering\\small\n",
         "\\begin{tabular}{lrr}\n\\toprule\n快速路径 & 1t & 4t \\\\\\n\\midrule\n",
         paste(frows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n\\end{minipage}\n")
}

batch_table <- function() {
  rows <- paste0(batch$n, " & ", fmt(batch$loop_ns / batch$n), " & ", fmt(batch$batch_ns / batch$n),
                 " & ", batch$loop_bytes, " & ", batch$batch_bytes, " \\\\")
  paste0("\\begin{tabular}{rrrrr}\\n\\toprule\\nN & 逐元素 ns/点 & 批量 ns/点 & 逐元素 B/点 & 批量 B/点 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

reuse_table <- function() {
  d <- batchscale[batchscale$section == "reuse", ]
  rows <- paste0(d$a, " & ", d$b, " & ", fmt(d$ns), " & ", sprintf("%.0f", d$bytes),
                 " & ", sprintf("%.1f", d$extra), " \\\\")
  paste0("\\begin{tabular}{rlrrr}\\n\\toprule\\nN & 变体 & ns/点 & B/点 & GB/s \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

joml_table <- function() {
  sel <- c("E3,rotate", "JOML,rotateAlloc", "P3,applyPoint",
           "JOML,matrixTransformAlloc", "E3,toQuaternion", "JOML,toMatrix4f")
  rows <- vapply(sel, function(s) { p <- strsplit(s, ",", fixed = TRUE)[[1]]
    v <- joml$ns[joml$path == p[1] & joml$op == p[2]][1]
    sprintf("\\texttt{%s} & \\texttt{%s} & %s \\\\", tex(p[1]), tex(p[2]), fmt(v)) }, "")
  paste0("\\begin{tabular}{llr}\\n\\toprule\\n路径 & 运算 & ns/op \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

render_table <- function() {
  rows <- paste0("\\texttt{", tex(render$variant), "} & ", render$n, " & ",
                 fmt(render$ns_per_point), " & ", fmt(render$p99_ns / 1e6), " & ",
                 sprintf("%.0f", render$bytes / 1e6), " \\\\")
  paste0("\\begin{tabular}{lrrrr}\\n\\toprule\\n路径 & 点数 & ns/点 & p99 帧 ms & 分配 MB/帧 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

robot_table <- function() {
  rows <- paste0("\\texttt{", tex(robot$op), "} & ", fmt(robot$ns), " \\\\")
  paste0("\\begin{tabular}{lr}\\n\\toprule\\n机器人求解 & ns/op \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

multicore_table <- function() {
  d <- multicore[multicore$threads %in% c(1, 4, 8, 14, 20), ]
  rows <- paste0("\\texttt{", tex(d$workload), "} & ", d$threads, " & ",
                 sprintf("%.2f", d$speedup), " & ", sprintf("%.0f\\%%", d$efficiency * 100), " \\\\")
  paste0("\\begin{tabular}{lrrr}\\n\\toprule\\n负载 & 线程 & 加速比 & 效率 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

fuzz_table <- function() {
  rows <- character()
  ds <- fuzz[fuzz$cat == "denseSparse" & fuzz$a %in% c(4, 8, 10), ]
  for (n in sort(unique(ds$a))) {
    rows <- c(rows, sprintf("\\texttt{denseSparse} & n=%s, gp/sandwich & %.2e \\\\",
                            n, max(ds$value[ds$a == n])))
  }
  od <- fuzz[fuzz$cat == "order", ]
  for (i in seq_len(nrow(od))) rows <- c(rows,
    sprintf("\\texttt{order} & n=%s & %.2e \\\\", od$a[i], od$value[i]))
  cd <- fuzz[fuzz$cat == "cond" & fuzz$b %in% c("1e-03", "1e-08", "1e-12"), ]
  for (i in seq_len(nrow(cd))) rows <- c(rows,
    sprintf("\\texttt{cond} & eps=%s & %.2e \\\\", cd$b[i], cd$value[i]))
  sm <- fuzz[fuzz$cat == "sum", ]
  if (nrow(sm)) rows <- c(rows, sprintf("\\texttt{sum} & double vs BigDecimal & %.2e \\\\", sm$value[1]))
  jm <- fuzz[fuzz$cat == "joml", ]
  if (nrow(jm)) rows <- c(rows, sprintf("\\texttt{joml} & float vs double & %.2e \\\\", jm$value[1]))
  paste0("\\begin{tabular}{llr}\\n\\toprule\\n类别 & 条件 & 结果 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

gabm_table <- function() {
  if (is.null(gabm)) return("")
  rows <- paste0("\\texttt{", tex(gabm$metric), "} & ", gabm$value, " \\\\")
  paste0("\\begin{tabular}{lr}\\n\\toprule\\nga-benchmark 指标 & 值 \\\\\\\\\n\\midrule\\n",
         paste(rows, collapse = "\n"), "\\n\\bottomrule\\n\\end{tabular}\\n")
}

template <- paste(readLines("template2.tex", warn = FALSE), collapse = "\n")
fill <- function(text, key, value) {
  parts <- strsplit(text, key, fixed = TRUE)[[1]]
  if (length(parts) == 1) stop("missing placeholder: ", key)
  paste(parts, collapse = value)
}
vals <- list(
  "%%SUMMARY%%" = summary,
  "%%SWEEP%%" = fixnl(sweep_table()),
  "%%FAST%%" = fixnl(fast_table()),
  "%%SUPPORT%%" = fixnl(support_table()),
  "%%SUBMASK%%" = fixnl(submask_table()),
  "%%BIGSPARSE%%" = fixnl(bigsparse_table()),
  "%%THRESH%%" = fixnl(threshold_table()),
  "%%SIGTAB%%" = fixnl(sigtab_table()),
  "%%ROBUST%%" = fixnl(robust_table()),
  "%%STRESS%%" = fixnl(stress_table()),
  "%%JMH%%" = fixnl(jmh_table()),
  "%%BATCH%%" = fixnl(batch_table()),
  "%%REUSE%%" = fixnl(reuse_table()),
  "%%MULTICORE%%" = fixnl(multicore_table()),
  "%%JOML%%" = fixnl(joml_table()),
  "%%RENDER%%" = fixnl(render_table()),
  "%%ROBOT%%" = fixnl(robot_table()),
  "%%FUZZ%%" = fixnl(fuzz_table()),
  "%%GABM%%" = fixnl(gabm_table()),
  "%%GP10%%" = sprintf("%.2f", gp10 / 1e6),
  "%%PAIR10%%" = sprintf("%.2f", pair10),
  "%%P3GP%%" = fmt(support$ns[support$kind == "P3"]),
  "%%E3GP%%" = fmt(support$ns[support$kind == "E3"]),
  "%%SUPPORT10%%" = sprintf("%.1f", support$ns[support$kind == "sparse10" & support$k == 10] / 1000),
  "%%OP10%%" = sprintf("%.0f", op10 / 1000),
  "%%OP10OLD%%" = sprintf("%.0f", sub10$before / 1000),
  "%%OPGAIN%%" = sprintf("%.1f", sub10$before / sub10$after),
  "%%SPARSE10%%" = sprintf("%.2f", sp10 / 1000),
  "%%SPARSEGAIN%%" = sprintf("%.0f", gp10 / sp10),
  "%%SPARSE14MS%%" = sprintf("%.0f", sp14$ns / 1e6),
  "%%SAND14%%" = sprintf("%.0f", sand14$ns / 1e6),
  "%%ALLOC1M%%" = sprintf("%.1f", alloc1m$ns),
  "%%REUSE1M%%" = sprintf("%.1f", reuse1m$ns),
  "%%DATE%%" = meta_val("date"),
  "%%CPU%%" = meta_val("cpu"),
  "%%JDK%%" = meta_val("jdk"),
  "%%JOMLVERSION%%" = meta_val("joml"),
  "%%TARGET%%" = meta_val("target"),
  "%%WARMUP%%" = meta_val("warmup"),
  "%%ROUNDS%%" = meta_val("rounds")
)
for (k in names(vals)) template <- fill(template, k, vals[[k]])
writeLines(template, "performance2.tex", useBytes = TRUE)
unlink(c("performance2.pdf", "performance2.aux"))
for (pass in 1:2) {
  log <- system2("xelatex", c("-interaction=nonstopmode", "-halt-on-error", "performance2.tex"),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  if (!is.null(status) && status != 0) { cat(log, sep = "\n"); stop("xelatex failed") }
}
if (!file.exists("performance2.pdf")) stop("xelatex produced no PDF")
file.copy("performance2.pdf", "../PERFORMANCE.pdf", overwrite = TRUE)
cat("wrote ../PERFORMANCE.pdf (2.0)\n")

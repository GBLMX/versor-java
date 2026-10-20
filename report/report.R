#!/usr/bin/env Rscript
# Builds PERFORMANCE.pdf from the CSVs in report/data.
# Requires: R (base + cairo), xelatex, Source Han Serif/Sans CN, DejaVu Sans Mono.
# Run: cd report && Rscript report.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
    getwd()
}
setwd(script_dir)
dir.create("figures", showWarnings = FALSE)

sweep <- read.csv("data/sweep.csv", stringsAsFactors = FALSE)
batch <- read.csv("data/batch.csv", stringsAsFactors = FALSE)
stress <- read.csv("data/stress.csv", stringsAsFactors = FALSE)
joml <- read.csv("data/joml.csv", stringsAsFactors = FALSE)
jit <- read.csv("data/jit.csv", stringsAsFactors = FALSE)
sparse <- read.csv("data/sparse.csv", stringsAsFactors = FALSE)
cliff <- read.csv("data/cliff.csv", stringsAsFactors = FALSE)
bigsparse <- read.csv("data/big-sparse.csv", stringsAsFactors = FALSE)
jmh <- read.csv("data/jmh.csv", stringsAsFactors = FALSE)
render <- read.csv("data/render.csv", stringsAsFactors = FALSE)
scale <- read.csv("data/scale.csv", stringsAsFactors = FALSE)
robot <- read.csv("data/robot.csv", stringsAsFactors = FALSE)
robust <- read.csv("data/robust.csv", stringsAsFactors = FALSE)
threshold <- read.csv("data/threshold.csv", stringsAsFactors = FALSE)
sigtab <- read.csv("data/table.csv", stringsAsFactors = FALSE)
batchscale <- read.csv("data/batchscale.csv", stringsAsFactors = FALSE)
fuzz <- read.csv("data/fuzz.csv", stringsAsFactors = FALSE)
meta <- read.csv("data/meta.csv", stringsAsFactors = FALSE)
meta_val <- function(key) meta$value[meta$key == key]

CJK <- "Source Han Sans CN"
fmt <- function(x) {
    ifelse(x >= 1000, formatC(round(x), format = "d", big.mark = ","), sprintf("%.1f", x))
}
tex <- function(x) gsub("_", "\\\\_", x)

# ------------------------------------------------------------------ figures

cairo_pdf("figures/fig-gp.pdf", width = 6.2, height = 4.0, family = CJK)
gp <- sweep[sweep$op == "gp", ]
gp <- gp[order(gp$n), ]
plot(gp$n, gp$ns, log = "xy", type = "b", pch = 19, lwd = 2, col = "#1f5fa8",
     xlab = "代数维度 n", ylab = "单次几何积 (ns, 对数轴)",
     main = "稠密几何积随维度按 4^n 增长")
ref <- gp$ns[1] * 4^(gp$n - gp$n[1])
lines(gp$n, ref, lty = 2, lwd = 1.5, col = "grey40")
legend("topleft", c("实测 gp", "参考 4^n"), lty = c(1, 2), pch = c(19, NA),
       col = c("#1f5fa8", "grey40"), bty = "n")
dev.off()

cairo_pdf("figures/fig-perpair.pdf", width = 6.2, height = 3.6, family = CJK)
per <- gp$ns / 4^gp$n
plot(gp$n, per, type = "b", pch = 19, lwd = 2, col = "#b45309", ylim = c(0, 1.3),
     xlab = "代数维度 n", ylab = "每对 blade 的耗时 (ns)",
     main = "单对成本近似恒定：瓶颈是 4^n 对循环")
abline(h = mean(per), lty = 2, col = "grey40")
dev.off()

cairo_pdf("figures/fig-batch.pdf", width = 6.2, height = 4.0, family = CJK)
plot(batch$n, batch$loop_ns, log = "xy", type = "b", pch = 19, lwd = 2,
     col = "#1f5fa8", xlab = "批量大小 N", ylab = "每点耗时 (ns, 对数轴)",
     main = "逐元素与批量变换（P3 motor）")
lines(batch$n, batch$batch_ns, type = "b", pch = 17, lwd = 2, col = "#b45309")
legend("bottomright", c("Multivector[] 逐元素", "MultivectorBatch (SoA)"),
       pch = c(19, 17), col = c("#1f5fa8", "#b45309"), bty = "n")
dev.off()

cairo_pdf("figures/fig-joml.pdf", width = 6.2, height = 3.8, family = CJK)
sel <- c("E3,rotate", "JOML,rotateAlloc", "P3,applyPoint",
         "JOML,matrixTransformAlloc", "E3,toQuaternion", "JOML,toMatrix4f")
labels <- c("E3 rotate", "JOML quat", "P3 motor×point",
            "JOML matrix", "E3→quat", "motor→M4")
vals <- sapply(sel, function(s) {
    p <- strsplit(s, ",", fixed = TRUE)[[1]]
    joml$ns[joml$path == p[1] & joml$op == p[2]][1]
})
par(mar = c(4, 10, 3, 1))
barplot(vals, names.arg = labels, horiz = TRUE, las = 1, log = "x",
        col = "#3f7f5f", xlab = "ns/op（对数轴）",
        main = "GA 与 JOML：运算成本与转换成本")
dev.off()

cairo_pdf("figures/fig-jit.pdf", width = 6.2, height = 3.6, family = CJK)
m <- rbind(jit$c2_ns, jit$c1_ns)
colnames(m) <- c("gp", "inverse", "dense sand", "P3 sand")
par(mar = c(5, 5, 3, 1))
barplot(m, beside = TRUE, log = "y", col = c("#1f5fa8", "#c2410c"),
        ylab = "ns/op（对数轴）", main = "C2 与仅 C1 的 JIT 对比",
        las = 1, cex.names = 0.95)
legend("topright", c("C2", "仅 C1"), fill = c("#1f5fa8", "#c2410c"), bty = "n")
dev.off()

cairo_pdf("figures/fig-sparse.pdf", width = 6.2, height = 4.0, family = CJK)
sp <- bigsparse[bigsparse$op == "gp" & bigsparse$k > 0, ]
plot(NA, xlim = range(sp$k), ylim = range(sp$ns), log = "xy",
     xlab = "非零 blade 数 k (对数轴)", ylab = "单次 gp (ns, 对数轴)",
     main = "高维稀疏几何积：成本随 k 而非 4^n")
cols <- c("#1f5fa8", "#b45309", "#3f7f5f")
for (i in seq_along(c(10, 12, 14))) {
    d <- sp[sp$n == c(10, 12, 14)[i], ]
    d <- d[order(d$k), ]
    lines(d$k, d$ns, type = "b", pch = 19, col = cols[i], lwd = 2)
}
dense10 <- bigsparse$ns[bigsparse$op == "gp" & bigsparse$n == 10 & bigsparse$k == -1]
abline(h = dense10, lty = 2, col = "grey40")
legend("topleft", c(paste0("n=", c(10, 12, 14)), "n=10 稠密"),
       col = c(cols, "grey40"), lty = c(1, 1, 1, 2), pch = c(19, 19, 19, NA), bty = "n")
dev.off()

cairo_pdf("figures/fig-render.pdf", width = 6.2, height = 3.6, family = CJK)
par(mar = c(4, 8, 3, 1))
barplot(rev(render$ns_per_point), names.arg = rev(render$variant), horiz = TRUE,
        las = 1, log = "x", col = "#3f7f5f", xlab = "每点耗时 (ns, 对数轴)",
        main = "每帧 20 万点的变换路径")
dev.off()

cairo_pdf("figures/fig-scale.pdf", width = 6.2, height = 3.6, family = CJK)
plot(scale$threads, scale$points_per_sec / 1e6, type = "b", pch = 19, lwd = 2,
     col = "#1f5fa8", xlab = "线程数", ylab = "吞吐 (M 点/秒)",
     main = "批量 sandwich 的线程扩展（20 万点）")
dev.off()

cairo_pdf("figures/fig-threshold.pdf", width = 6.2, height = 3.8, family = CJK)
th <- threshold[threshold$n == 14, ]
plot(th$k, th$hash_ns, log = "xy", type = "b", pch = 19, col = "#b45309", lwd = 2,
     xlab = "非零数 k (对数轴)", ylab = "gp (ns, 对数轴)",
     main = "n=14 策略交叉：哈希 vs 稠密累加")
lines(th$k, th$dense_ns, type = "b", pch = 17, col = "#1f5fa8", lwd = 2)
lines(th$k, th$auto_ns, type = "l", lty = 2, col = "grey30", lwd = 1.5)
legend("topleft", c("哈希", "稠密累加", "自动选择"),
       pch = c(19, 17, NA), lty = c(1, 1, 2), col = c("#b45309", "#1f5fa8", "grey30"), bty = "n")
dev.off()

cairo_pdf("figures/fig-cond.pdf", width = 6.2, height = 3.6, family = CJK)
cd <- fuzz[fuzz$cat == "cond", ]
plot(as.numeric(cd$a), pmax(cd$value, 1e-30), log = "xy", type = "b", pch = 19,
     col = "#b45309", lwd = 2, xlab = "扰动 eps（对数轴）",
     ylab = "残差 |A x - 1|（对数轴）",
     main = "inverseGeneral 近奇异残差")
dev.off()

cairo_pdf("figures/fig-reuse.pdf", width = 6.2, height = 3.6, family = CJK)
r1 <- batchscale[batchscale$section == "reuse" & batchscale$a == "1000000", ]
par(mar = c(4, 7, 3, 1))
barplot(rev(r1$ns), names.arg = rev(r1$b), horiz = TRUE, las = 1,
        col = "#3f7f5f", xlab = "ns/点",
        main = "100 万点：分配 vs 输出复用")
dev.off()

# ------------------------------------------------------------------- tables

sweep_table <- function() {
    ops <- c("gp", "op", "ip", "rightContraction", "commutator", "reverse",
             "dual", "norm2", "unit", "inverse", "sandwich",
             "inverseGeneral", "exp")
    dims <- sort(unique(sweep$n))
    head <- paste0("\\textbf{运算} & ", paste0("\\textbf{n=", dims, "}", collapse = " & "), " \\\\")
    rows <- vapply(ops, function(op) {
        cells <- vapply(dims, function(n) {
            v <- sweep$ns[sweep$op == op & sweep$n == n]
            if (!length(v)) "--" else fmt(v[1])
        }, "")
        paste0("\\texttt{", tex(op), "} & ", paste(cells, collapse = " & "), " \\\\")
    }, "")
    paste0("\\resizebox{\\textwidth}{!}{%\n",
           "{\\scriptsize\\setlength{\\tabcolsep}{3pt}\n",
           "\\begin{tabular}{l", strrep("r", length(dims)), "}\n\\toprule\n",
           head, "\n\\midrule\n", paste(rows, collapse = "\n"),
           "\n\\bottomrule\n\\end{tabular}}}\n")
}

pair_table <- function() {
    d <- data.frame(n = gp$n, pairs = formatC(4^gp$n, format = "d", big.mark = ","),
                    ns = sprintf("%.3f", per))
    rows <- paste0(d$n, " & ", 2^d$n, " & ", d$pairs, " & ", d$ns, " \\\\")
    paste0("\\begin{tabular}{rrrr}\n\\toprule\n",
           "n & 系数 $2^n$ & blade 对 $4^n$ & ns/pair \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

batch_table <- function() {
    rows <- paste0(batch$n, " & ", fmt(batch$loop_ns), " & ", fmt(batch$batch_ns),
                   " & ", batch$loop_bytes, " & ", batch$batch_bytes, " \\\\")
    paste0("\\begin{tabular}{rrrrr}\n\\toprule\n",
           "N & 逐元素 ns/点 & 批量 ns/点 & 逐元素 B/点 & 批量 B/点 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

cliff_table <- function() {
    rows <- vapply(seq_len(nrow(cliff)), function(i) {
        r <- cliff[i, ]
        sprintf("%s & %d & %s & %s & %.1f$\\times$ \\\\", tex(r$op), r$n,
                fmt(r$before_ns), fmt(r$after_ns), r$before_ns / r$after_ns)
    }, "")
    paste0("\\begin{tabular}{lr rrr}\n\\toprule\n",
           "运算 & n & 修复前 ns & 修复后 ns & 提升 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

stress_table <- function() {
    rows <- vapply(seq_len(nrow(stress)), function(i) {
        r <- stress[i, ]
        m3 <- if (is.na(r$matrix3_inverse_ns) || r$matrix3_inverse_ns == "") "--" else fmt(r$matrix3_inverse_ns)
        sprintf("%d & %s & %s & %s \\\\", r$n, fmt(r$dense_sandwich_ns),
                fmt(r$inverse_general_ns), m3)
    }, "")
    paste0("\\begin{tabular}{rrrr}\n\\toprule\n",
           "n & 稠密 sandwich & inverseGeneral & 3×3 P3 矩阵求逆 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

joml_table <- function() {
    rows <- paste0("\\texttt{", tex(joml$path), "} & \\texttt{", tex(joml$op),
                   "} & ", sprintf("%.2f", joml$ns), " \\\\")
    paste0("\\begin{tabular}{llr}\n\\toprule\n路径 & 运算 & ns/op \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

jit_table <- function() {
    rows <- paste0("\\texttt{", tex(jit$op), "} & ", fmt(jit$c2_ns), " & ",
                   fmt(jit$c1_ns), " & ", sprintf("%.2f", jit$c1_ns / jit$c2_ns),
                   "$\\times$ \\\\")
    paste0("\\begin{tabular}{lrrr}\n\\toprule\n运算 & C2 ns/op & 仅 C1 ns/op & C2 加速 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

sparse_table <- function() {
    rows <- paste0("\\texttt{", tex(sparse$metric), "} & ", sprintf("%.2f", sparse$ns), " \\\\")
    paste0("\\begin{tabular}{lr}\n\\toprule\n指标 & ns/op \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

bigsparse_table <- function() {
    pick <- bigsparse[bigsparse$op == "gp" & bigsparse$n %in% c(10, 12, 14) & bigsparse$k %in% c(2, 8, 10, 32, 128, 512), ]
    rows <- vapply(seq_len(nrow(pick)), function(i) {
        r <- pick[i, ]
        sprintf("%d & %d & %s & %s & %s \\\\", r$n, r$k,
                fmt(bigsparse$ns[bigsparse$n == r$n & bigsparse$k == r$k & bigsparse$op == "gp"]),
                fmt(bigsparse$ns[bigsparse$n == r$n & bigsparse$k == r$k & bigsparse$op == "sandwich"]),
                sprintf("%.2f", r$ns_per_pair))
    }, "")
    dens <- bigsparse[bigsparse$k == -1 & bigsparse$n == 10, ]
    drows <- paste0("10 & 稠密 & ", fmt(dens$ns[dens$op == "gp"]), " & ",
                    fmt(dens$ns[dens$op == "sandwich"]), " & -- \\\\")
    paste0("\\begin{tabular}{rrrrr}\n\\toprule\n",
           "n & k & gp ns & sandwich ns & ns/pair \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\midrule\n",
           paste(drows, collapse = "\n"),
           "\n\\bottomrule\n\\end{tabular}\n")
}

jmh_table <- function() {
    d <- jmh[jmh$bench == "dense", ]
    dims <- sort(unique(d$n))
    rows <- vapply(dims, function(n) {
        v <- function(op) sprintf("%.1f", d$ns[d$n == n & d$op == op][1])
        sprintf("%d & %s & %s & %s \\\\", n, v("gp"), v("sandwich"), v("inverse"))
    }, "")
    f <- jmh[jmh$bench == "fast" & jmh$threads == 1, ]
    frows <- paste0("\\texttt{", tex(f$op), "} & ", sprintf("%.1f", f$ns), " \\\\")
    paste0("\\begin{tabular}{rrrr}\n\\toprule\n",
           "n & gp & sandwich & inverse \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"),
           "\n\\bottomrule\n\\end{tabular}\n\\vspace{0.5em}\n",
           "\\begin{tabular}{lr}\n\\toprule\n快速路径 / 批量 & ns/op \\\\\n\\midrule\n",
           paste(frows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

render_table <- function() {
    rows <- vapply(seq_len(nrow(render)), function(i) {
        r <- render[i, ]
        sprintf("\\texttt{%s} & %d & %s & %s & %.0f \\\\", tex(r$variant), r$n,
                fmt(r$ns_per_point), fmt(r$p99_ns / 1e6), r$bytes)
    }, "")
    paste0("\\begin{tabular}{lrrrr}\n\\toprule\n",
           "路径 & 点数 & ns/点 & p99 帧 (ms) & 分配 B/帧 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

scale_table <- function() {
    base <- scale$points_per_sec[scale$threads == 1]
    rows <- paste0(scale$threads, " & ", sprintf("%.1f", scale$points_per_sec / 1e6),
                   " & ", sprintf("%.2f", scale$ns_per_point),
                   " & ", sprintf("%.2f$\\times$", scale$points_per_sec / base), " \\\\")
    paste0("\\begin{tabular}{rrrr}\n\\toprule\n",
           "线程 & M 点/秒 & ns/点 & 相对 1 线程 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

robot_table <- function() {
    rows <- paste0("\\texttt{", tex(robot$op), "} & ", sprintf("%.2f", robot$ns), " \\\\")
    paste0("\\begin{tabular}{lr}\n\\toprule\n机器人求解路径 & ns/op \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

robust_table <- function() {
    pick <- robust[robust$k %in% c(8, 32) & robust$kind %in% c("uniform", "versor", "pgaPoint")
        & robust$n %in% c(10, 12, 14), ]
    rows <- vapply(seq_len(nrow(pick)), function(i) {
        r <- pick[i, ]
        sprintf("%d & \\texttt{%s} & %d & \\texttt{%s} & %s & %.0f \\\\", r$n, r$kind,
                r$k, r$op, fmt(r$ns), r$p99)
    }, "")
    paste0("\\begin{tabular}{rlrlrr}\n\\toprule\n",
           "n & 结构 & k & 运算 & ns & p99 ns \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

threshold_table <- function() {
    rows <- paste0(threshold$n, " & ", threshold$k, " & ", fmt(threshold$hash_ns), " & ",
                   fmt(threshold$dense_ns), " & ", fmt(threshold$auto_ns), " & ",
                   threshold$choice, " & ", threshold$correct, " \\\\")
    paste0("\\begin{tabular}{rrrrrll}\n\\toprule\n",
           "n & k & hash & dense & auto & 选择 & 判定 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

sigtab_table <- function() {
    rows <- paste0(sigtab$n, " & ", sigtab$threads, " & ", sprintf("%.1f", sigtab$first_us),
                   " & ", sprintf("%.1f", sigtab$rebuild_us), " & ",
                   ifelse(sigtab$bytes == 0, "--", sprintf("%.1f MiB", sigtab$bytes / 1048576)),
                   " \\\\")
    paste0("\\begin{tabular}{rrrrl}\n\\toprule\n",
           "n & 线程 & 首次乘积 us & 重建 us & 表内存 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

batchscale_table <- function() {
    r1 <- batchscale[batchscale$section == "reuse", ]
    rows <- paste0(r1$a, " & ", r1$b, " & ", fmt(r1$ns), " & ",
                   sprintf("%.0f", r1$bytes), " & ", sprintf("%.1f", r1$extra), " \\\\")
    paste0("\\begin{tabular}{rlrrr}\n\\toprule\n",
           "N & 变体 & ns/点 & B/点 & GB/s \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

batchscale_fig_table <- function() {
    r1 <- batchscale[batchscale$section == "scale", ]
    rows <- paste0(r1$a, " & ", sprintf("%.1f", r1$extra / 1e6), " & ",
                   sprintf("%.2f", r1$ns), " \\\\")
    paste0("\\begin{tabular}{rrr}\n\\toprule\n",
           "线程 & M 点/秒 & ns/点 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

fuzz_table <- function() {
    rows <- paste0("\\texttt{", tex(fuzz$cat), "} & ", tex(fuzz$a), " & ",
                   sprintf("%.2e", fuzz$value), " \\\\")
    paste0("\\begin{tabular}{llr}\n\\toprule\n类别 & 指标 & 值 \\\\\n\\midrule\n",
           paste(rows, collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

# --------------------------------------------------------------------- tex

template <- paste(readLines("template.tex", warn = FALSE), collapse = "\n")
fill <- function(text, key, value) {
    parts <- strsplit(text, key, fixed = TRUE)[[1]]
    if (length(parts) == 1) stop("missing placeholder: ", key)
    paste(parts, collapse = value)
}

template <- fill(template, "%%SWEEP%%", sweep_table())
template <- fill(template, "%%BIGSPARSE%%", bigsparse_table())
template <- fill(template, "%%JMH%%", jmh_table())
template <- fill(template, "%%RENDER%%", render_table())
template <- fill(template, "%%SCALE%%", scale_table())
template <- fill(template, "%%ROBOT%%", robot_table())
template <- fill(template, "%%ROBUST%%", robust_table())
template <- fill(template, "%%THRESHOLD%%", threshold_table())
template <- fill(template, "%%SIGTAB%%", sigtab_table())
template <- fill(template, "%%BATCHSCALE%%", batchscale_table())
template <- fill(template, "%%BATCHSCALEFIG%%", batchscale_fig_table())
template <- fill(template, "%%FUZZ%%", fuzz_table())
template <- fill(template, "%%PAIRS%%", pair_table())
template <- fill(template, "%%CLIFF%%", cliff_table())
template <- fill(template, "%%BATCH%%", batch_table())
template <- fill(template, "%%SPARSE%%", sparse_table())
template <- fill(template, "%%JOML%%", joml_table())
template <- fill(template, "%%STRESS%%", stress_table())
template <- fill(template, "%%JIT%%", jit_table())
template <- fill(template, "%%CPU%%", meta_val("cpu"))
template <- fill(template, "%%JDK%%", meta_val("jdk"))
template <- fill(template, "%%JOMLVERSION%%", meta_val("joml"))
template <- fill(template, "%%DATE%%", meta_val("date"))
template <- fill(template, "%%TARGET%%", meta_val("target"))
template <- fill(template, "%%WARMUP%%", meta_val("warmup"))
template <- fill(template, "%%ROUNDS%%", meta_val("rounds"))
writeLines(template, "performance.tex", useBytes = TRUE)

for (pass in 1:2) {
    log <- system2("xelatex",
                   c("-interaction=nonstopmode", "-halt-on-error", "performance.tex"),
                   stdout = TRUE, stderr = TRUE)
    if (!file.exists("performance.pdf")) {
        cat(log, sep = "\n")
        stop("xelatex failed")
    }
}
file.copy("performance.pdf", "../PERFORMANCE.pdf", overwrite = TRUE)
cat("wrote ../PERFORMANCE.pdf\n")

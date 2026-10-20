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

# --------------------------------------------------------------------- tex

template <- paste(readLines("template.tex", warn = FALSE), collapse = "\n")
fill <- function(text, key, value) {
    parts <- strsplit(text, key, fixed = TRUE)[[1]]
    if (length(parts) == 1) stop("missing placeholder: ", key)
    paste(parts, collapse = value)
}

template <- fill(template, "%%SWEEP%%", sweep_table())
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

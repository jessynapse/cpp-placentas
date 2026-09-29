# ==============================================================================
# 34-figure1-flow.R
# CPP Placental Pathology → Infantile Hemangioma
# September 2026
#
# Figure 1: participant flow diagram, from outputs/flow_counts.csv (33).
# Output: manuscript/figure1_flow.png (600 dpi) and .pdf
# ==============================================================================

library(tidyverse)

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()
fc <- read_csv(file.path(root, "outputs", "flow_counts.csv"), show_col_types = FALSE)
n <- setNames(fc$n, fc$step); f <- function(x) format(x, big.mark = ",")

main <- c(
  sprintf("Pregnancies (infants) in the\nCollaborative Perinatal Project\nn = %s", f(n[1])),
  sprintf("Infants with a placental\npathology record\nn = %s", f(n[3])),
  sprintf("Analytic sample: singleton pregnancies\nwith MVM and AI data\nn = %s", f(n[7])),
  sprintf("Infants with an observed\nIH outcome at one year\nn = %s\n(%s definite IH cases)", f(n[10]), f(n[11])))
side <- c(
  sprintf("Excluded: no placental\npathology record (n = %s)", f(n[2])),
  sprintf("Excluded (n = %s)\n  Multiple gestation (n = %s)\n  Missing MVM score (n = %s)\n  Missing AI (n = %s)",
          f(n[4] + n[5] + n[6]), f(n[4]), f(n[5]), f(n[6])),
  sprintf("Excluded (n = %s)\n  Died before one year (n = %s)\n  No one-year IH examination (n = %s)",
          f(n[8] + n[9]), f(n[8]), f(n[9])))

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10))
  ys <- c(9, 6.4, 3.8, 1.1); w <- 2.1; h <- 0.85
  for (i in 1:4) {
    rect(3.2 - w, ys[i] - h, 3.2 + w, ys[i] + h, lwd = 1.2)
    text(3.2, ys[i], main[i], cex = 0.78)
    if (i < 4) arrows(3.2, ys[i] - h, 3.2, ys[i + 1] + h, length = 0.08, lwd = 1.2)
  }
  for (i in 1:3) {
    ym <- (ys[i] + ys[i + 1]) / 2
    segments(3.2, ym, 5.8, ym, lwd = 1.2)
    rect(5.8, ym - 0.75, 9.9, ym + 0.75, lwd = 1.2)
    text(5.95, ym, side[i], adj = 0, cex = 0.72)
  }
}
png(file.path(root, "manuscript", "figure1_flow.png"), width = 6.5, height = 6, units = "in", res = 600)
draw(); dev.off()
pdf(file.path(root, "manuscript", "figure1_flow.pdf"), width = 6.5, height = 6)
draw(); dev.off()
cat("Saved: manuscript/figure1_flow.png, .pdf\n")

# ==============================================================================
# 44-dag-grouped.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Grouped traditional DAG (supersedes 43). Covariates are grouped into
# domain nodes (each node is a vector of variables, all adjusted, so
# boxed). Arrows within a domain are documented in the eTable of DAG
# references rather than drawn.
#   Sociodemographic Factors -> Maternal Health and Behaviour, Study Site
#   every boxed node -> Placental MVM or AI and Infantile Haemangioma
#   Study Site -> selection (site-specific follow-up)
#   Placental MVM or AI -> mediators -> Infantile Haemangioma
#   Placental MVM or AI and mediators -> selection
# Labels in title case. Output: manuscript/efigure_dag.png (600 dpi), .pdf
# ==============================================================================

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()

fam <- tryCatch(
  if (any(grepl("^Calibri$", systemfonts::system_fonts()$family))) "Calibri" else "sans",
  error = function(e) "sans")
dgrey <- "grey30"

# nodes: x, y, half-width, half-height
N <- list(
  soc  = c(2.9, 8.3, 2.25, 0.95), hlth = c(8.6, 8.3, 2.5, 0.95),
  site = c(14.3, 8.3, 1.3, 0.45), sex  = c(8.6, 5.1, 1.0, 0.33),
  A    = c(2.5, 3.6, 1.7, 0.45),   Y    = c(14.4, 3.6, 1.7, 0.45),
  M    = c(8.4, 1.95, 2.6, 0.5),   S    = c(14.4, 0.2, 1.7, 0.45)
)
head <- c(soc = "Sociodemographic Factors", hlth = "Maternal Health and Behaviour",
          site = "Study Site", sex = "Infant Sex")
body <- c(soc = "Race/Ethnicity, Maternal Age,\nEducation, Family Income,\nMarital Status",
          hlth = "Pre-Pregnancy BMI, Pre-Pregnancy Diabetes,\nChronic Hypertension, Parity,\nSmoking During Pregnancy")

narr <- function(f, t, bend = 0, col = "black", lty = 1, lwd = 1.1, pad = 0.05, len = 0.07) {
  x0 <- f[1]; y0 <- f[2]; x1 <- t[1]; y1 <- t[2]
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  mx <- (x0 + x1) / 2 - bend * (y1 - y0) / d
  my <- (y0 + y1) / 2 + bend * (x1 - x0) / d
  u <- seq(0, 1, length.out = 500)
  bx <- (1 - u)^2 * x0 + 2 * (1 - u) * u * mx + u^2 * x1
  by <- (1 - u)^2 * y0 + 2 * (1 - u) * u * my + u^2 * y1
  ins <- function(n) abs(bx - n[1]) <= n[3] + pad & abs(by - n[2]) <= n[4] + pad
  k <- !ins(f) & !ins(t); bx <- bx[k]; by <- by[k]; n <- length(bx)
  lines(bx, by, col = col, lty = lty, lwd = lwd)
  arrows(bx[n - 1], by[n - 1], bx[n], by[n], length = len, col = col, lty = lty, lwd = lwd)
}
a <- function(f, t, ...) narr(N[[f]], N[[t]], ...)

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = fam)
  plot.new(); plot.window(xlim = c(0.4, 16.9), ylim = c(-0.45, 10.0))

  # among confounder domains
  a("soc", "hlth"); a("soc", "site", bend = 2.8)
  # confounders -> exposure and outcome
  for (k in c("soc", "hlth", "site", "sex")) { a(k, "A"); a(k, "Y") }
  # effect of interest, mediators
  a("A", "Y", lwd = 2.8, len = 0.11)
  a("A", "M"); a("M", "Y")
  # selection
  a("A", "S", bend = -3.0, col = dgrey); a("M", "S", col = dgrey)
  a("site", "S", bend = 4.0, col = dgrey)
  # unmeasured
  U <- c(8.6, 4.2, 0.3, 0.3)
  narr(U, N$A, lty = 2, bend = 0.2); narr(U, N$Y, lty = 2, bend = -0.2)

  # nodes
  for (k in c("soc", "hlth")) {
    n <- N[[k]]
    rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], col = "white", lwd = 1.2)
    text(n[1], n[2] + 0.55, head[[k]], cex = 0.78, font = 2)
    text(n[1], n[2] - 0.2, body[[k]], cex = 0.68)
  }
  for (k in c("site", "sex")) {
    n <- N[[k]]
    rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], col = "white", lwd = 1.2)
    text(n[1], n[2], head[[k]], cex = 0.78, font = 2)
  }
  n <- N$S
  rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], col = "white", lwd = 1.8)
  text(n[1], n[2], "Survived and Examined\nat 12 Months (S = 1)", cex = 0.72)
  n <- N$M
  rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], col = "white", border = NA)
  text(n[1], n[2], "Pre-Eclampsia, Gestational Age at Delivery,\nFoetal Growth (Mediators, Not Adjusted)", cex = 0.7)
  text(N$A[1], N$A[2], "Placental MVM or AI", cex = 0.95, font = 2)
  text(N$Y[1], N$Y[2], "Infantile Haemangioma", cex = 0.95, font = 2)
  symbols(U[1], U[2], circles = 0.27, inches = FALSE, add = TRUE, lty = 2, bg = "white")
  text(U[1], U[2], "U", cex = 0.8)
  text(U[1] + 0.4, U[2] + 0.28, "Unmeasured Common Causes", adj = 0, cex = 0.62, col = dgrey)

  legend(0.45, 1.2, bty = "n", cex = 0.62, y.intersp = 1.3,
         legend = c("Boxed: conditioned on (adjusted for, or restricted to S = 1)",
                    "Dashed: unmeasured",
                    "Grey: causes of selection, addressed with IPCW"),
         lty = c(NA, 2, 1), col = c("black", "black", dgrey), pch = c(0, NA, NA))
}

png(file.path(root, "manuscript", "efigure_dag.png"), width = 11, height = 6.6,
    units = "in", res = 600, type = "cairo")
draw(); dev.off()
cairo_pdf(file.path(root, "manuscript", "efigure_dag.pdf"), width = 11, height = 6.6)
draw(); dev.off()
cat("Saved: manuscript/efigure_dag.png, .pdf (font:", fam, ")\n")

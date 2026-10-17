# ==============================================================================
# 43-dag-traditional.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Traditional single-panel DAG for Paediatric and Perinatal Epidemiology.
#   - exposure (left) -> outcome (right), effect of interest in bold
#   - each adjusted covariate is its own node, boxed (box = conditioned on)
#   - covariate -> exposure and covariate -> outcome arrows in light grey
#   - arrows among covariates in black (each supported by the literature,
#     see eTable of DAG references)
#   - mediators (not adjusted) below the main arrow
#   - unmeasured common causes (U) dashed
#   - selection (survived and examined at 12 months) boxed, with its causes
# Output: manuscript/efigure_dag.png (600 dpi) and .pdf
# ==============================================================================

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()

fam <- tryCatch(
  if (any(grepl("^Calibri$", systemfonts::system_fonts()$family))) "Calibri" else "sans",
  error = function(e) "sans")
lgrey <- "grey62"; dgrey <- "grey30"

# nodes: x, y, half-width, half-height
N <- list(
  race  = c(2.0, 9.3, 1.25, 0.3), educ = c(5.3, 9.3, 1.0, 0.3),
  marit = c(8.5, 9.3, 1.1, 0.3),  inc  = c(11.7, 9.3, 1.1, 0.3),
  site  = c(14.9, 9.3, 0.95, 0.3),
  age   = c(2.0, 7.5, 1.05, 0.3), par  = c(4.75, 7.5, 0.7, 0.3),
  bmi   = c(7.35, 7.5, 1.2, 0.3), dm   = c(10.15, 7.5, 1.25, 0.3),
  htn   = c(13.0, 7.5, 1.25, 0.3), smoke = c(15.35, 7.5, 0.75, 0.3),
  sex   = c(8.5, 5.55, 0.85, 0.3),
  A     = c(2.3, 3.6, 1.55, 0.45), Y = c(14.7, 3.6, 1.35, 0.45),
  M     = c(7.6, 1.95, 2.1, 0.45),
  U     = c(8.5, 3.95, 0.0, 0.0),
  S     = c(14.6, 0.2, 1.5, 0.42)
)
lab <- c(race = "Race and ethnicity", educ = "Education", marit = "Marital status",
         inc = "Family income", site = "Study site", age = "Maternal age",
         par = "Parity", bmi = "Pre-pregnancy BMI", dm = "Pre-pregnancy diabetes",
         htn = "Chronic hypertension", smoke = "Smoking", sex = "Infant sex",
         A = "Placental MVM or AI", Y = "Infantile haemangioma",
         M = "Pre-eclampsia, gestational age at delivery,\nfoetal growth (mediators, not adjusted)",
         S = "Survived and examined\nat 12 months (S = 1)")
covs <- c("race", "educ", "marit", "inc", "site", "age", "par", "bmi", "dm", "htn", "smoke", "sex")

narr <- function(f, t, bend = 0, col = "black", lty = 1, lwd = 1, pad = 0.05, len = 0.06) {
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
  plot.new(); plot.window(xlim = c(0.4, 17.4), ylim = c(-0.45, 9.95))

  # covariates -> exposure and outcome (light grey)
  for (k in covs) { a(k, "A", col = lgrey, lwd = 0.7); a(k, "Y", col = lgrey, lwd = 0.7) }

  # among covariates (black)
  a("race", "educ"); a("marit", "inc"); a("inc", "site")
  a("race", "marit", bend = 0.5); a("race", "inc", bend = 0.75)
  a("race", "site", bend = 1.0); a("educ", "inc", bend = 0.45)
  a("educ", "age"); a("marit", "age")
  a("age", "par"); a("par", "bmi")
  a("race", "bmi"); a("educ", "bmi"); a("inc", "bmi"); a("inc", "smoke")
  a("bmi", "dm"); a("bmi", "htn", bend = -1.0)
  a("age", "dm", bend = -1.7); a("age", "htn", bend = -2.2)

  # exposure, mediators, outcome, selection, U
  a("A", "Y", lwd = 2.6, len = 0.1)
  a("A", "M", lwd = 1.1); a("M", "Y", lwd = 1.1)
  a("A", "S", bend = -2.6, col = dgrey); a("M", "S", col = dgrey)
  a("site", "S", bend = 5.0, col = dgrey)
  narr(c(8.5, 4.55, 0.3, 0.3), N$A, lty = 2, bend = 0.25)
  narr(c(8.5, 4.55, 0.3, 0.3), N$Y, lty = 2, bend = -0.25)

  # nodes
  for (k in names(N)) {
    n <- N[[k]]; if (k == "U") next
    boxed <- k %in% c(covs, "S")
    rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], col = "white",
         border = if (boxed) "black" else NA, lwd = if (k == "S") 1.8 else 1.1)
    text(n[1], n[2], lab[[k]], cex = if (k %in% c("A", "Y")) 0.9 else 0.68,
         font = if (k %in% c("A", "Y")) 2 else 1)
  }
  symbols(8.5, 4.55, circles = 0.27, inches = FALSE, add = TRUE, lty = 2, bg = "white")
  text(8.5, 4.55, "U", cex = 0.8)
  text(9.0, 4.55, "Unmeasured common causes", adj = 0, cex = 0.6, col = dgrey)

  legend(0.45, 1.25, bty = "n", cex = 0.6, y.intersp = 1.25,
         legend = c("Boxed: conditioned on (adjusted for, or restricted to S = 1)",
                    "Grey: every boxed covariate affects exposure and outcome",
                    "Dashed: unmeasured"),
         lty = c(NA, 1, 2), col = c("black", lgrey, "black"), pch = c(0, NA, NA))
}

png(file.path(root, "manuscript", "efigure_dag.png"), width = 11, height = 7.2,
    units = "in", res = 600, type = "cairo")
draw(); dev.off()
cairo_pdf(file.path(root, "manuscript", "efigure_dag.pdf"), width = 11, height = 7.2)
draw(); dev.off()
cat("Saved: manuscript/efigure_dag.png, .pdf (font:", fam, ")\n")

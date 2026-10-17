# ==============================================================================
# 42-dag-faithful.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Causal DAG drawn to the conventions of Hernan and Robins (What If):
#   - a box around a node means the analysis conditions on it (adjustment
#     or restriction). No grouping boxes.
#   - an omitted arrow is a strong assumption, so arrows between covariate
#     tiers are saturated (any variable may affect any later variable)
#   - unmeasured common causes are drawn (U)
#   - true and measured exposure and outcome are separate nodes (A, A*,
#     Y, Y*), following Hernan and Cole (2009), because study site acts
#     mainly through how pathology and IH were recorded
#
# Panel A: temporal tiers of the baseline covariates L. Because every node
#   in L is adjusted for, arrows among them do not change the adjustment
#   set, and they are shown for completeness.
# Panel B: structure used for identification.
#   L  baseline covariates (Panel A), boxed (adjusted)
#   V  study site, boxed (adjusted). Caused by L (residence). Affects true
#      A and Y (local population, environment, care) and the recording of
#      both (pathologist reading practice -> A*, paediatric ascertainment
#      -> Y*). V -> A* and V -> Y* is a dependent-error path that only
#      adjustment for V blocks.
#   U  unmeasured common causes of A and Y (residual confounding, E-value)
#   A  true placental MVM or AI.  A* recorded MVM or AI (the exposure used)
#   M  pre-eclampsia, gestational age, foetal growth (descendants of A,
#      not adjusted, so the estimate is a total effect)
#   Y  true IH.  Y* definite IH recorded at 12 months (the outcome used)
#   S  survived to one year and examined at 12 months, boxed (restricted).
#      Caused by L, V, A and M. Because M -> S and M -> Y, restricting to
#      S = 1 opens A -> S <- M -> Y. IPCW with A, L and V only does not
#      block this path. Censoring weights that also include gestational
#      age and birthweight would.
# Output: manuscript/efigure_dag.png (600 dpi) and .pdf
# ==============================================================================

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()

fam <- tryCatch(
  if (any(grepl("^Calibri$", systemfonts::system_fonts()$family))) "Calibri" else "sans",
  error = function(e) "sans")
red <- "#8B0000"; grey <- "grey35"

# node: text with an optional conditioning box
node <- function(x, y, lab, sub = NULL, boxed = FALSE, dashed = FALSE, w = 0.75, h = 0.38, cex = 0.95) {
  if (boxed) rect(x - w, y - h, x + w, y + h, lwd = 1.6)
  if (dashed) symbols(x, y, circles = 0.42, inches = FALSE, add = TRUE, lty = 2, lwd = 1.1)
  text(x, y, lab, cex = cex, font = 2)
  if (!is.null(sub)) text(x, y - h - 0.28, sub, cex = 0.58, col = grey)
}
# straight or curved arrow between points, trimmed at both ends
carr <- function(x0, y0, x1, y1, bend = 0, trim0 = 0.45, trim1 = 0.45,
                 col = "black", lty = 1, lwd = 1.1) {
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  mx <- (x0 + x1) / 2 - bend * (y1 - y0) / d
  my <- (y0 + y1) / 2 + bend * (x1 - x0) / d
  t <- seq(0, 1, length.out = 60)
  bx <- (1 - t)^2 * x0 + 2 * (1 - t) * t * mx + t^2 * x1
  by <- (1 - t)^2 * y0 + 2 * (1 - t) * t * my + t^2 * y1
  arc <- c(0, cumsum(sqrt(diff(bx)^2 + diff(by)^2)))
  keep <- arc >= trim0 & arc <= max(arc) - trim1
  bx <- bx[keep]; by <- by[keep]; n <- length(bx)
  lines(bx, by, col = col, lty = lty, lwd = lwd)
  arrows(bx[n - 1], by[n - 1], bx[n], by[n], length = 0.07, col = col, lty = lty, lwd = lwd)
}

panel_a <- function() {
  par(mar = c(0.2, 0.2, 1.4, 0.2), family = fam)
  plot.new(); plot.window(xlim = c(0, 14.6), ylim = c(0, 6.2))
  title("A. Baseline covariates (L), in temporal order", adj = 0, cex.main = 0.95, font.main = 2)
  tx <- c(1.4, 4.9, 9.2, 13.2)
  hd <- c("Fixed", "Social position", "Before this pregnancy", "During pregnancy")
  text(tx, 5.85, hd, cex = 0.7, font = 3, col = grey)
  bw <- 1.3; bh <- 0.3
  b <- function(x, y, lab) { rect(x - bw, y - bh, x + bw, y + bh, lwd = 1.3); text(x, y, lab, cex = 0.68) }
  b(tx[1], 3.6, "Race and ethnicity")
  b(tx[2], 5.0, "Education"); b(tx[2], 3.85, "Marital status")
  b(tx[2], 2.7, "Family income"); b(tx[2], 1.55, "Maternal age")
  b(tx[3], 5.0, "Parity"); b(tx[3], 3.85, "Pre-pregnancy BMI")
  b(tx[3], 2.7, "Pre-pregnancy diabetes"); b(tx[3], 1.55, "Chronic hypertension")
  b(tx[4], 3.6, "Smoking")
  b(tx[1], 0.45, "Infant sex"); text(tx[1] + bw + 0.15, 0.45,
    "Determined at conception: no causes in L. Affects A and Y.", adj = 0, cex = 0.6, font = 3, col = grey)
  # saturated arrows between tiers
  for (i in 1:3) {
    xa <- tx[i] + bw + c(0.12, 0.85, 0.75)[i]; xb <- tx[i + 1] - bw - 0.12
    polygon(c(xa, xb - 0.25, xb - 0.25, xb, xb - 0.25, xb - 0.25, xa),
            c(3.45, 3.45, 3.3, 3.6, 3.9, 3.75, 3.75), col = "grey80", border = "grey50")
  }
  text(7.3, 0.95, "Grey arrows: every variable may affect every variable in a later tier",
       cex = 0.6, font = 3, col = grey)
  # arrows within tiers (arcs bend outward from the column)
  carr(tx[2] - bw, 5.0, tx[2] - bw, 2.7, bend = -0.45, trim0 = 0.02, trim1 = 0.02)   # education -> income
  carr(tx[2], 3.85 - bh, tx[2], 2.7 + bh, trim0 = 0.02, trim1 = 0.02)                # marital -> income
  carr(tx[2] + bw, 5.0, tx[2] + bw, 1.55, bend = 0.7, trim0 = 0.02, trim1 = 0.02)    # education -> age
  carr(tx[2] + bw, 3.85, tx[2] + bw, 1.55, bend = 0.35, trim0 = 0.02, trim1 = 0.02)  # marital -> age
  carr(tx[3] - bw, 5.0, tx[3] - bw, 3.85, bend = -0.3, trim0 = 0.02, trim1 = 0.02)   # parity -> BMI
  carr(tx[3], 3.85 - bh, tx[3], 2.7 + bh, trim0 = 0.02, trim1 = 0.02)                # BMI -> diabetes
  carr(tx[3] + bw, 3.85, tx[3] + bw, 1.55, bend = 0.45, trim0 = 0.02, trim1 = 0.02)  # BMI -> chronic HTN
}

# arrow between two nodes (x, y, half-width, half-height), clipped at both borders
narr <- function(f, t, bend = 0, col = "black", lty = 1, lwd = 1.1, pad = 0.06) {
  x0 <- f[1]; y0 <- f[2]; x1 <- t[1]; y1 <- t[2]
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  mx <- (x0 + x1) / 2 - bend * (y1 - y0) / d
  my <- (y0 + y1) / 2 + bend * (x1 - x0) / d
  u <- seq(0, 1, length.out = 400)
  bx <- (1 - u)^2 * x0 + 2 * (1 - u) * u * mx + u^2 * x1
  by <- (1 - u)^2 * y0 + 2 * (1 - u) * u * my + u^2 * y1
  inside <- function(n) abs(bx - n[1]) <= n[3] + pad & abs(by - n[2]) <= n[4] + pad
  keep <- !inside(f) & !inside(t)
  bx <- bx[keep]; by <- by[keep]; n <- length(bx)
  lines(bx, by, col = col, lty = lty, lwd = lwd)
  arrows(bx[n - 1], by[n - 1], bx[n], by[n], length = 0.07, col = col, lty = lty, lwd = lwd)
}

panel_b <- function() {
  par(mar = c(0.2, 0.2, 1.4, 0.2), family = fam)
  plot.new(); plot.window(xlim = c(0, 13.4), ylim = c(-2.3, 7.6))
  title("B. Causal structure for the effect of placental pathology on IH", adj = 0,
        cex.main = 0.95, font.main = 2)
  # x, y, half-width, half-height
  P <- list(U  = c(8.6, 6.9, 1.75, 0.42), L = c(1.35, 3.6, 1.25, 0.42),
            V  = c(1.35, 0.6, 1.0, 0.32),  A = c(5.2, 3.6, 1.05, 0.3),
            As = c(5.0, 0.6, 1.25, 0.42), M = c(9.6, 2.2, 1.55, 0.42),
            Y  = c(11.9, 3.6, 0.8, 0.3),   Ys = c(12.0, 0.6, 1.3, 0.42),
            S  = c(8.0, 0.6, 1.2, 0.42))
  a <- function(f, t, ...) narr(P[[f]], P[[t]], ...)
  # confounding and measurement
  a("L", "A"); a("L", "Y", bend = 1.25); a("L", "V")
  a("V", "A"); a("V", "As"); a("V", "Y", bend = 1.0)
  a("V", "Ys", bend = -2.3)
  a("U", "A", lty = 2); a("U", "Y", lty = 2)
  a("A", "As"); a("Y", "Ys")
  # effect of interest and mediators
  a("A", "Y", col = red, lwd = 2.4); a("A", "M"); a("M", "Y")
  # selection
  a("M", "S", col = grey); a("A", "S", col = grey); a("L", "S", col = grey)
  a("V", "S", bend = -1.35, col = grey)

  nd <- function(k, lab, boxed = FALSE, dashed = FALSE, cex = 0.66) {
    n <- P[[k]]
    if (boxed) rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], lwd = 1.6, col = "white")
    if (dashed) rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], lty = 2, col = "white")
    if (!boxed && !dashed) rect(n[1] - n[3], n[2] - n[4], n[1] + n[3], n[2] + n[4], col = "white", border = NA)
    text(n[1], n[2], lab, cex = cex)
  }
  nd("U", "U: unmeasured common causes\n(e.g. diet, infection, environment)", dashed = TRUE)
  nd("L", "L: baseline covariates\n(Panel A)", boxed = TRUE)
  nd("V", "V: study site", boxed = TRUE)
  nd("A", "A: true MVM or AI", cex = 0.72)
  nd("As", "A*: recorded MVM or AI\n(exposure analysed)")
  nd("M", "M: pre-eclampsia, gestational\nage, foetal growth (not adjusted)")
  nd("Y", "Y: true IH", cex = 0.72)
  nd("Ys", "Y*: definite IH recorded\nat 12 months (outcome analysed)")
  nd("S", "S = 1: survived and\nexamined at 12 months", boxed = TRUE)
  text(8.55, 3.83, "Effect of interest (total effect)", cex = 0.62, font = 3, col = red)
  text(0.05, -2.2, paste(
    "Boxed: conditioned on (adjusted or restricted). Dashed: unmeasured. Grey arrows: causes of selection.",
    "V -> A* and V -> Y* (reading practice, ascertainment) is blocked only by adjusting for V.",
    sep = "  "), adj = 0, cex = 0.56, col = grey)
}

draw <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(1, 1.25))
  panel_a(); panel_b()
}

png(file.path(root, "manuscript", "efigure_dag.png"), width = 10, height = 9.6,
    units = "in", res = 600, type = "cairo")
draw(); dev.off()
cairo_pdf(file.path(root, "manuscript", "efigure_dag.pdf"), width = 10, height = 9.6)
draw(); dev.off()
cat("Saved: manuscript/efigure_dag.png, .pdf (font:", fam, ")\n")

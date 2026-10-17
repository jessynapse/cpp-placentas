# ==============================================================================
# 41-dag-detailed.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Detailed DAG (supersedes the layouts in 39 and 40). Every adjusted
# covariate is its own node, laid out in tiers of assumed causal order:
#   tier 1 (background): study site, race and ethnicity, maternal age
#   tier 2 (social):     education, marital status, income
#   tier 3 (behaviour and reproductive history): smoking, parity
#   tier 4 (pre-pregnancy health): BMI, diabetes, chronic hypertension
#   infant sex: no maternal causes
# Assumed links among covariates:
#   site -> race (catchment), education, income
#   race -> education, income, BMI, chronic hypertension
#   age -> marital status, parity, diabetes, chronic hypertension
#   education -> income, smoking, BMI
#   marital status -> income, smoking
#   income -> smoking
#   parity -> BMI
#   BMI -> diabetes, chronic hypertension
# All covariates -> placental pathology and IH (the confounder box).
# Site also acts through pathologist reading practice (exposure) and
# paediatric ascertainment (outcome).
# Not adjusted (dashed): gestational age and foetal growth, downstream of
# placental pathology. They also affect survival (S).
# S = 1: survived to one year and examined at 12 months, affected by the
# covariates, placental pathology and gestational age. Addressed with
# stabilised IPCW (denominator: MVM, AI, all 11 covariates and site).
# Output: manuscript/efigure_dag.png (600 dpi) and .pdf
# ==============================================================================

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()

fam <- tryCatch(
  if (any(grepl("^Calibri$", systemfonts::system_fonts()$family))) "Calibri" else "sans",
  error = function(e) "sans")

blue_fill <- "#DAE8FC"; blue_line <- "#6C8EBF"
orng_fill <- "#FFE6CC"; orng_line <- "#D79B00"
red <- "#8B0000"; grey <- "grey30"

# ---- layout grid ----
cx <- c(1.75, 5.25, 8.75, 12.25)          # tier columns, half-width 1.3
ry <- c(9.5, 8.0, 6.5, 5.0)               # rows, half-height 0.42
W <- 1.3; H <- 0.42
L <- function(i) cx[i] - W; R <- function(i) cx[i] + W
Tp <- function(j) ry[j] + H; Bt <- function(j) ry[j] - H

nodes <- list(
  site  = c(1, 1, "Study site"),
  race  = c(1, 2, "Race and ethnicity"),
  age   = c(1, 3, "Maternal age"),
  sex   = c(1, 4, "Infant sex"),
  educ  = c(2, 1, "Education"),
  inc   = c(2, 2, "Family income"),
  marit = c(2, 3, "Marital status"),
  smoke = c(3, 2, "Smoking during\npregnancy"),
  par   = c(3, 4, "Parity"),
  bmi   = c(4, 1, "Pre-pregnancy BMI"),
  dm    = c(4, 2, "Pre-pregnancy\ndiabetes"),
  htn   = c(4, 3, "Chronic\nhypertension")
)

rbox <- function(x, y, w, h, lab, fill = blue_fill, border = blue_line,
                 cex = 0.72, lty = 1, lwd = 1.1, r = 0.15) {
  t <- seq(0, pi / 2, length.out = 12)
  qx <- c(x + w - r, x - w + r, x - w + r, x + w - r)
  qy <- c(y + h - r, y + h - r, y - h + r, y - h + r)
  px <- c(qx[1] + r * cos(t), qx[2] + r * cos(t + pi / 2),
          qx[3] + r * cos(t + pi), qx[4] + r * cos(t + 3 * pi / 2))
  py <- c(qy[1] + r * sin(t), qy[2] + r * sin(t + pi / 2),
          qy[3] + r * sin(t + pi), qy[4] + r * sin(t + 3 * pi / 2))
  polygon(px, py, col = fill, border = border, lty = lty, lwd = lwd)
  if (!is.null(lab)) text(x, y, lab, cex = cex)
}
arr <- function(x0, y0, x1, y1, col = "black", lty = 1, lwd = 1.1)
  arrows(x0, y0, x1, y1, length = 0.08, col = col, lty = lty, lwd = lwd)
# right-angled path through waypoints, arrowhead on the last segment
pth <- function(x, y, col = "grey15", lty = 1, lwd = 0.95) {
  n <- length(x)
  if (n > 2) lines(x[-n], y[-n], col = col, lty = lty, lwd = lwd)
  arrows(x[n - 1], y[n - 1], x[n], y[n], length = 0.06, col = col, lty = lty, lwd = lwd)
}

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = fam)
  plot.new(); plot.window(xlim = c(0, 18.2), ylim = c(-0.35, 10.95))

  # ---- confounder box ----
  bx <- c(0.15, 14.35, 4.4, 10.7)          # left, right, bottom, top
  rect(bx[1], bx[3], bx[2], bx[4], lwd = 1)
  text(bx[1] + 0.15, 10.47, "Confounders (all adjusted in the primary model)",
       adj = 0, cex = 0.72, col = "#1F3A93")
  tiers <- c("Background", "Social", "Behaviour and reproductive history", "Pre-pregnancy health")
  text(cx, 10.17, tiers, cex = 0.58, font = 3, col = grey)

  # links among confounders
  pth(c(R(1), L(2)), c(9.65, 9.65))                                # site -> education
  pth(c(cx[1], cx[1]), c(Bt(1), Tp(2)))                            # site -> race
  pth(c(R(1), 3.75, 3.75, L(2)), c(9.25, 9.25, 8.15, 8.15))       # site -> income
  lines(c(R(1), 3.35, 3.35), c(7.9, 7.9, 9.45), col = "grey15", lwd = 0.95)  # race trunk
  lines(c(3.35, 3.35), c(7.9, 7.25), col = "grey15", lwd = 0.95)
  pth(c(3.35, L(2)), c(9.45, 9.45))                                # race -> education
  pth(c(3.35, L(2)), c(7.9, 7.9))                                  # race -> income
  pth(c(3.35, 11.45, 11.45), c(8.8, 8.8, Bt(1)))                   # race -> BMI
  pth(c(3.35, 10.4, 10.4, L(4)), c(7.25, 7.25, 6.6, 6.6))          # race -> chronic HTN
  pth(c(R(1), L(2)), c(6.5, 6.5))                                  # age -> marital status
  lines(c(3.6, 3.6, 13.95, 13.95), c(6.5, 5.75, 5.75, 8.0), col = "grey15", lwd = 0.95)  # age trunk
  pth(c(cx[3], cx[3]), c(5.75, Tp(4)))                             # age -> parity
  pth(c(cx[4] + 0.4, cx[4] + 0.4), c(5.75, Bt(3)))                 # age -> chronic HTN
  pth(c(13.95, R(4)), c(8.0, 8.0))                                 # age -> diabetes
  pth(c(cx[2], cx[2]), c(Bt(1), Tp(2)))                            # education -> income
  pth(c(R(2), L(4)), c(9.6, 9.6))                                  # education -> BMI
  pth(c(R(2), 7.0, 7.0, L(3)), c(9.35, 9.35, 8.15, 8.15))          # education -> smoking
  pth(c(R(2), L(3)), c(8.0, 8.0))                                  # income -> smoking
  pth(c(cx[2], cx[2]), c(Tp(3), Bt(2)))                            # marital status -> income
  pth(c(R(2), 7.0, 7.0, L(3)), c(6.5, 6.5, 7.85, 7.85))            # marital status -> smoking
  pth(c(R(3), 10.75, 10.75, L(4)), c(5.0, 5.0, 9.38, 9.38))        # parity -> BMI
  pth(c(cx[4], cx[4]), c(Bt(1), Tp(2)))                            # BMI -> diabetes
  pth(c(R(4), 14.15, 14.15, R(4)), c(9.5, 9.5, 6.5, 6.5))          # BMI -> chronic HTN

  for (k in names(nodes)) {
    n <- nodes[[k]]; i <- as.numeric(n[1]); j <- as.numeric(n[2])
    rbox(cx[i], ry[j], W, H, n[3])
  }
  text(R(1) + 0.1, 4.75, "(no maternal causes)", adj = 0, cex = 0.56, font = 3, col = grey)

  # ---- exposure, outcome, mediators ----
  ex <- c(3.1, 2.75, 2.75, 0.95); ou <- c(11.4, 2.75, 2.75, 0.95)
  rbox(ex[1], ex[2], ex[3], ex[4], NULL, fill = orng_fill, border = orng_line, lwd = 2)
  text(ex[1], ex[2] + 0.6, "Exposure: placental pathology", font = 2, cex = 0.78)
  text(ex[1], ex[2] - 0.22,
       "(1) MVM, any lesion (and score, 0 to 8)\n(2) AI, any maternal or foetal\nneutrophilic inflammation\n(3) Joint: MVM only, AI only, both",
       cex = 0.66)
  rbox(ou[1], ou[2], ou[3], ou[4], NULL, fill = orng_fill, border = orng_line, lwd = 2)
  text(ou[1], ou[2] + 0.6, "Outcome: infantile haemangioma", font = 2, cex = 0.78)
  text(ou[1], ou[2] - 0.15,
       "Definite IH recorded at the\n12-month paediatric examination\n(suspect IH counted as absent)",
       cex = 0.66)
  me <- c(7.25, 0.35, 2.45, 0.48)
  rbox(me[1], me[2], me[3], me[4],
       "Gestational age at delivery and foetal growth\n(potential mediators, not adjusted)",
       fill = "white", border = "grey35", cex = 0.64, lty = 2)

  arr(ex[1], bx[3], ex[1], ex[2] + ex[4])                      # confounders -> exposure
  arr(ou[1], bx[3], ou[1], ou[2] + ou[4])                      # confounders -> outcome
  text(ex[1] + 0.12, 4.02, "incl. site: pathologist\nreading practice", adj = 0, cex = 0.56, font = 3, col = grey)
  text(ou[1] + 0.12, 4.02, "incl. site: paediatric\nascertainment of IH", adj = 0, cex = 0.56, font = 3, col = grey)
  arr(ex[1] + ex[3], ex[2], ou[1] - ou[3], ou[2], col = red, lwd = 2.2)  # exposure -> outcome
  text(7.25, 2.98, "Association of interest", cex = 0.62, font = 3, col = red)
  text(7.25, 2.53, "(total effect)", cex = 0.58, font = 3, col = red)
  arr(ex[1] + 1.6, ex[2] - ex[4], me[1] - 1.2, me[2] + me[4], lty = 2)  # exposure -> mediators
  arr(me[1] + 1.2, me[2] + me[4], ou[1] - 1.6, ou[2] - ou[4], lty = 2)  # mediators -> outcome

  # ---- selection node ----
  sn <- c(16.35, 7.55, 1.7, 1.35)
  rbox(sn[1], sn[2], sn[3], sn[4], NULL, fill = "white", border = "black", lwd = 2.2)
  text(sn[1], sn[2] + 0.85, "S = 1: survived to one year\nand examined at 12 months", font = 2, cex = 0.64)
  text(sn[1], sn[2] - 0.35,
       "40,700 of 45,268 included\n(1,834 died, 2,734 not examined)\n\nAddressed with stabilised IPCW\n(MVM, AI, 11 covariates, site)",
       cex = 0.58)
  arr(bx[2], sn[2], sn[1] - sn[3], sn[2])                      # confounders -> S
  lines(c(ex[1] - 2.0, ex[1] - 2.0, 17.65), c(ex[2] - ex[4], -0.25, -0.25), lwd = 1.1)  # exposure -> S
  arr(17.65, -0.25, 17.65, sn[2] - sn[4])
  lines(c(me[1] + me[3], 16.95), c(me[2] - 0.12, me[2] - 0.12), lty = 2, lwd = 1.1)  # mediators -> S
  arr(16.95, me[2] - 0.12, 16.95, sn[2] - sn[4], lty = 2)
  arr(15.4, sn[2] - sn[4], ou[1] + ou[3], ou[2] + 0.45, lty = 3)  # S -> IH observed
  text(14.75, 3.95, "IH observed\nonly if S = 1", adj = 0, cex = 0.56, font = 3)
}

png(file.path(root, "manuscript", "efigure_dag.png"), width = 11.5, height = 7.4,
    units = "in", res = 600, type = "cairo")
draw(); dev.off()
cairo_pdf(file.path(root, "manuscript", "efigure_dag.pdf"), width = 11.5, height = 7.4)
draw(); dev.off()
cat("Saved: manuscript/efigure_dag.png, .pdf (font:", fam, ")\n")

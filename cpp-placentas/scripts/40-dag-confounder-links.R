# ==============================================================================
# 40-dag-confounder-links.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Redrawn DAG (replaces the layout from 39-dag.R) that shows the assumed
# links among the confounders, in the style of the ECHO Aim 1 DAG:
#   - confounder box (all adjusted in the primary model), with arrows
#     between covariates:
#       study site -> race and ethnicity, socioeconomic position
#         (site determines the catchment population)
#       race and ethnicity -> socioeconomic position, cardiometabolic health
#       socioeconomic position -> smoking, cardiometabolic health
#       maternal age -> parity, cardiometabolic health
#       parity -> cardiometabolic health
#       infant sex has no maternal causes
#   - confounder box -> placental pathology and IH
#   - placental pathology -> IH (association of interest)
#   - gestational age and foetal growth are potential mediators (dashed,
#     not adjusted)
#   - S (survived and examined at one year) is the selection node,
#     affected by the confounders and placental pathology, and addressed
#     with inverse probability of censoring weights
# Output: manuscript/efigure_dag.png (600 dpi) and .pdf
# ==============================================================================

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()

# Calibri for figure labels (journal style) when installed, otherwise sans
fam <- tryCatch(
  if (any(grepl("^Calibri$", systemfonts::system_fonts()$family))) "Calibri" else "sans",
  error = function(e) "sans")

blue_fill <- "#DAE8FC"; blue_line <- "#6C8EBF"
orng_fill <- "#FFE6CC"; orng_line <- "#D79B00"
red       <- "#8B0000"

# rounded box centred at (x, y) with half-width w and half-height h
rbox <- function(x, y, w, h, lab, fill = blue_fill, border = blue_line,
                 cex = 0.8, lty = 1, lwd = 1.2, font = 1, r = 0.18) {
  t <- seq(0, pi / 2, length.out = 12)
  cx <- c(x + w - r, x - w + r, x - w + r, x + w - r)
  cy <- c(y + h - r, y + h - r, y - h + r, y - h + r)
  px <- c(cx[1] + r * cos(t), cx[2] + r * cos(t + pi / 2),
          cx[3] + r * cos(t + pi), cx[4] + r * cos(t + 3 * pi / 2))
  py <- c(cy[1] + r * sin(t), cy[2] + r * sin(t + pi / 2),
          cy[3] + r * sin(t + pi), cy[4] + r * sin(t + 3 * pi / 2))
  polygon(px, py, col = fill, border = border, lty = lty, lwd = lwd)
  text(x, y, lab, cex = cex, font = font)
}
arr <- function(x0, y0, x1, y1, col = "black", lty = 1, lwd = 1.1)
  arrows(x0, y0, x1, y1, length = 0.08, col = col, lty = lty, lwd = lwd)
seg <- function(x, y, col = "black", lty = 1, lwd = 1.1) lines(x, y, col = col, lty = lty, lwd = lwd)

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = fam)
  plot.new(); plot.window(xlim = c(0, 13.6), ylim = c(0, 10.6))

  # ---- confounder box ----
  rect(0.2, 4.4, 9.8, 10.5, border = "black", lwd = 1)
  text(0.35, 10.28, "Confounders (all adjusted in the primary model)",
       adj = 0, cex = 0.72, col = "#1F3A93")

  # columns x = 2.1, 5.7, 8.4 | rows y = 9.2, 7.4, 5.6
  w1 <- 1.55; w3 <- 1.0; h <- 0.55
  rbox(2.1, 9.2, w1, h, "Study site")
  rbox(5.7, 9.2, w1, h, "Socioeconomic position\n(education, income,\nmarital status)", cex = 0.72)
  rbox(8.65, 9.2, w3, h, "Smoking")
  rbox(2.1, 7.4, w1, h, "Race and ethnicity")
  rbox(5.7, 7.4, w1, h, "Pre-pregnancy BMI,\ndiabetes, chronic\nhypertension", cex = 0.72)
  rbox(8.65, 7.4, w3, h, "Infant sex")
  rbox(2.1, 5.6, w1, h, "Maternal age")
  rbox(5.7, 5.6, w1, h, "Parity")

  # links among confounders
  arr(2.1 + w1, 9.2, 5.7 - w1, 9.2)               # site -> SES
  arr(2.1, 9.2 - h, 2.1, 7.4 + h)                 # site -> race
  arr(3.0, 7.4 + h, 4.6, 9.2 - h)                 # race -> SES
  arr(2.1 + w1, 7.4, 5.7 - w1, 7.4)               # race -> cardiometabolic
  arr(5.7 + w1, 9.2, 8.65 - w3, 9.2)               # SES -> smoking
  arr(5.7, 9.2 - h, 5.7, 7.4 + h)                 # SES -> cardiometabolic
  arr(2.1 + w1, 5.6, 5.7 - w1, 5.6)               # age -> parity
  arr(3.0, 5.6 + h, 4.6, 7.4 - h)                 # age -> cardiometabolic
  arr(5.7, 5.6 + h, 5.7, 7.4 - h)                 # parity -> cardiometabolic

  # ---- exposure, outcome, mediators, selection ----
  rbox(2.6, 2.9, 1.9, 0.75, "\nPlacental MVM, AI, or both",
       fill = orng_fill, border = orng_line, cex = 0.85, lwd = 2)
  text(2.6, 3.3, "Exposure", font = 2, cex = 0.85)
  rbox(8.4, 2.9, 1.9, 0.75, "\nInfantile haemangioma\nat one year",
       fill = orng_fill, border = orng_line, cex = 0.85, lwd = 2)
  text(8.4, 3.3, "Outcome", font = 2, cex = 0.85)
  rbox(5.5, 0.95, 1.75, 0.55, "Gestational age and foetal growth\n(potential mediators, not adjusted)",
       fill = "white", border = "grey35", cex = 0.68, lty = 2)
  rbox(11.4, 7.4, 1.35, 0.85, "S: survived and\nexamined at one year\n(addressed with\nIPCW)",
       fill = "white", border = "black", cex = 0.68, lwd = 2.2)

  arr(2.6, 4.4, 2.6, 2.9 + 0.75)                  # confounders -> exposure
  arr(8.4, 4.4, 8.4, 2.9 + 0.75)                  # confounders -> outcome
  arr(2.6 + 1.9, 2.9, 8.4 - 1.9, 2.9, col = red, lwd = 2.2)   # exposure -> outcome
  text(5.5, 3.12, "Association of interest", cex = 0.62, font = 3, col = red)
  arr(3.4, 2.9 - 0.75, 4.2, 0.95 + 0.55, lty = 2)  # exposure -> mediators
  arr(6.8, 0.95 + 0.55, 7.6, 2.9 - 0.75, lty = 2)  # mediators -> outcome
  arr(9.8, 7.4, 11.4 - 1.35, 7.4)                 # confounders -> S
  seg(c(1.2, 1.2, 12.3), c(2.9 - 0.75, 0.2, 0.2)) # exposure -> S (routed below)
  arr(12.3, 0.2, 12.3, 7.4 - 0.85)
  arr(11.4, 7.4 - 0.85, 8.4 + 1.9, 2.9 + 0.3, lty = 3)  # S -> IH observed only if S = 1
  text(12.45, 4.2, "IH observed\nonly if S = 1", cex = 0.6, font = 3, adj = 0)
}

png(file.path(root, "manuscript", "efigure_dag.png"), width = 9, height = 7.4,
    units = "in", res = 600, type = "cairo")
draw(); dev.off()
cairo_pdf(file.path(root, "manuscript", "efigure_dag.pdf"), width = 9, height = 7.4)
draw(); dev.off()
cat("Saved: manuscript/efigure_dag.png, .pdf (font:", fam, ")\n")

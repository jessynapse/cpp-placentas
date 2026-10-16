# ==============================================================================
# 39-dag.R
# CPP Placental Pathology -> Infantile Hemangioma
# October 2026
#
# Directed acyclic graph (DAG), recommended by Paediatric and Perinatal
# Epidemiology, showing the assumed causal structure:
#   - confounders (11 covariates) affect placental pathology and IH
#   - study site affects recorded placental diagnoses (pathologist reading
#     practice) and IH (population, paediatric ascertainment)
#   - gestational age and foetal growth are potential mediators (not adjusted)
#   - survival to and examination at one year (S) is the selection node,
#     addressed with inverse probability of censoring weights
# Output: manuscript/efigure_dag.png (600 dpi) and .pdf
# ==============================================================================

root <- normalizePath(file.path(getwd(), ".."))
if (basename(getwd()) != "scripts") root <- getwd()

box <- function(x, y, w, h, lab, cex = 0.72, lty = 1, lwd = 1.2, border = "black") {
  rect(x - w, y - h, x + w, y + h, lty = lty, lwd = lwd, border = border)
  text(x, y, lab, cex = cex)
}
arr <- function(x0, y0, x1, y1, lty = 1) arrows(x0, y0, x1, y1, length = 0.08, lwd = 1.1, lty = lty)

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10))
  # nodes
  box(1.7, 7.6, 1.55, 1.3, paste("Confounders", "Maternal age, BMI,", "race and ethnicity,",
      "education, income,", "marital status, smoking,", "parity, diabetes,", "chronic hypertension,",
      "infant sex", sep = "\n"), cex = 0.62)
  box(5, 9.2, 1.3, 0.55, "Study site\n(population, pathologist,\nascertainment)", cex = 0.62)
  box(3.4, 4.2, 1.25, 0.75, "Placental pathology\nMVM, AI, or both")
  box(8.4, 4.2, 1.1, 0.6, "Infantile\nhaemangioma")
  box(5.9, 1.6, 1.35, 0.6, "Gestational age and\nfoetal growth (mediators,\nnot adjusted)", cex = 0.6, lty = 2)
  box(8.4, 7.2, 1.1, 0.65, "S: survived and\nexamined at one year\n(IPW)", cex = 0.6, lwd = 2.2)
  # arrows
  arr(4.65, 4.2, 7.3, 4.2)                     # exposure -> IH
  arr(2.2, 6.3, 2.9, 4.95); arr(3.25, 7.2, 7.4, 4.75)   # confounders -> exposure, IH
  arr(4.4, 8.65, 3.7, 4.95); arr(5.6, 8.65, 8.1, 4.8)   # site -> exposure, IH
  arr(3.9, 3.45, 5.0, 2.2, lty = 2); arr(6.9, 2.2, 8.0, 3.6, lty = 2)  # mediation
  arr(4.4, 4.95, 7.4, 6.85); arr(3.25, 8.4, 7.3, 7.45)  # exposure, confounders -> S
  arr(8.4, 6.55, 8.4, 4.8, lty = 3)                     # S: IH observed only if S = 1
  text(5.9, 4.45, "Association of interest", cex = 0.58, font = 3)
}
png(file.path(root, "manuscript", "efigure_dag.png"), width = 7, height = 6, units = "in", res = 600)
draw(); dev.off()
pdf(file.path(root, "manuscript", "efigure_dag.pdf"), width = 7, height = 6)
draw(); dev.off()
cat("Saved: manuscript/efigure_dag.png, .pdf\n")

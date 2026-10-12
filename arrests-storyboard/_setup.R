# Shared data preparation, helpers and models for the Arrests exploration and storyboard.
# Packages: install.packages(c("flexdashboard", "carData", "ggplot2", "vcd", "dplyr", "tidyr"))

library(carData)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vcd)

data(Arrests)
A <- Arrests %>%
  mutate(rel = as.integer(released == "Yes"),
         agegrp = cut(age, c(11, 17, 20, 25, 30, 40, 70),
                      labels = c("12-17", "18-20", "21-25", "26-30", "31-40", "41-66")),
         period = ifelse(year <= 1999, "1997 to 1999", "2000 to 2002"))

theme_set(theme_minimal(base_size = 13) +
            theme(panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))
col_race <- c(Black = "#3B6FB6", White = "#E1A33B")
col_rel  <- c(No = "#C8553D", Yes = "#5B8E7D")

pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")

# Proportion released with a Wilson 95% CI, by any grouping
rel_by <- function(data, ...) {
  data %>%
    group_by(...) %>%
    summarise(n = n(), k = sum(rel), .groups = "drop") %>%
    mutate(p = k / n,
           z = 1.96,
           lo = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
           hi = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n))
}

overall <- mean(A$rel)
rate <- function(cond) mean(A$rel[cond])

# Models
fit_crude <- glm(released ~ colour, family = binomial, data = A)
fit_adj   <- glm(released ~ colour + year + age + sex + employed + citizen + checks,
                 family = binomial, data = A)
fit_int   <- update(fit_adj, . ~ . + colour:year + colour:age)
fit_main3 <- glm(released ~ colour + employed + citizen, family = binomial, data = A)
fit_cec   <- glm(released ~ colour * employed * citizen, family = binomial, data = A)

or_ci <- function(fit, term) {
  est <- coef(fit)[term]
  se  <- sqrt(vcov(fit)[term, term])
  exp(c(est, est - 1.96 * se, est + 1.96 * se))
}
or_crude <- or_ci(fit_crude, "colourWhite")
or_adj   <- or_ci(fit_adj, "colourWhite")
lrt_int  <- anova(fit_adj, fit_int, test = "Chisq")
lrt_cec  <- anova(fit_main3, fit_cec, test = "Chisq")
fmt_or   <- function(o) sprintf("%.2f (95%% CI %.2f to %.2f)", o[1], o[2], o[3])
fmt_p    <- function(p) ifelse(p < 1e-4, "p < 0.0001", sprintf("p = %.2f", p))

# ---- Plot functions (shared so the storyboard reuses exploration plots) ----

plot_release_by_colour <- function() {
  rel_by(A, colour) %>%
    ggplot(aes(x = colour, y = p, fill = colour)) +
    geom_col(width = 0.55, show.legend = FALSE) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12) +
    geom_text(aes(y = 0.05, label = paste0(pct(p, 1), "\nn = ", scales::comma(n))),
              vjust = 0, colour = "white", size = 5, fontface = "bold") +
    scale_fill_manual(values = col_race) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Percent released with a summons (95% CI)",
         title = "Release with a summons by colour")
}

plot_profile_by_colour <- function() {
  A %>%
    group_by(colour) %>%
    summarise(`Unemployed` = mean(employed == "No"),
              `Not a citizen` = mean(citizen == "No"),
              `3 or more police checks` = mean(checks >= 3),
              `No police checks` = mean(checks == 0),
              `Male` = mean(sex == "Male"),
              `Under 21` = mean(age < 21)) %>%
    pivot_longer(-colour, names_to = "trait", values_to = "p") %>%
    mutate(trait = factor(trait, levels = rev(c("Unemployed", "Not a citizen",
                                                "3 or more police checks", "No police checks",
                                                "Male", "Under 21")))) %>%
    ggplot(aes(x = p, y = trait)) +
    geom_line(aes(group = trait), colour = "grey75", linewidth = 1.3) +
    geom_point(aes(colour = colour), size = 4) +
    geom_text(aes(label = pct(p), colour = colour, vjust = ifelse(colour == "Black", -1, 2)),
              size = 3.6, show.legend = FALSE) +
    scale_colour_manual(values = col_race, name = "Colour") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Percent of arrestees with the characteristic", y = NULL,
         title = "Characteristics of Black and White arrestees")
}

plot_marginal_factors <- function() {
  marg <- bind_rows(lapply(c("colour", "sex", "employed", "citizen"), function(v) {
    rel_by(A, level = .data[[v]]) %>% mutate(variable = v, level = as.character(level))
  })) %>%
    mutate(variable = factor(variable, levels = c("colour", "employed", "citizen", "sex")),
           label = paste0(variable, ": ", level))
  ggplot(marg, aes(x = p, y = reorder(label, p))) +
    geom_vline(xintercept = overall, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = "grey40") +
    geom_point(aes(size = n), colour = "#2F4858") +
    geom_text(aes(label = paste0(pct(p, 1), "  (n = ", n, ")")), nudge_y = 0.33, size = 3.6) +
    facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(labels = scales::percent, limits = c(0.6, 0.95)) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = "Percent released with a summons (95% CI)", y = NULL,
         title = "Release by each binary characteristic",
         subtitle = "Dashed line is the overall release rate")
}

plot_rate_by <- function(d, xlab, ttl, ylim = c(0.3, 1)) {
  ggplot(d, aes(x = factor(x), y = p)) +
    geom_hline(yintercept = overall, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, colour = "grey40") +
    geom_point(aes(size = n), colour = "#2F4858") +
    geom_line(aes(group = 1), colour = "#2F4858", alpha = 0.5) +
    scale_y_continuous(labels = scales::percent, limits = ylim) +
    scale_size_area(max_size = 5, guide = "none") +
    labs(x = xlab, y = "% released", title = ttl)
}

plot_release_drivers <- function() {
  d <- bind_rows(
    rel_by(A, lev = employed) %>% mutate(var = "Employed", lev = as.character(lev)),
    rel_by(A, lev = citizen)  %>% mutate(var = "Citizen",  lev = as.character(lev)),
    A %>% mutate(lev = ifelse(checks >= 5, "5-6", as.character(checks))) %>%
      rel_by(lev) %>% mutate(var = "Number of police checks")
  ) %>%
    mutate(var = factor(var, levels = c("Employed", "Citizen", "Number of police checks")))
  ggplot(d, aes(x = lev, y = p)) +
    geom_hline(yintercept = overall, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, colour = "grey40") +
    geom_point(aes(size = n), colour = "#2F4858") +
    facet_grid(. ~ var, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(0.5, 1)) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(x = NULL, y = "Percent released (95% CI)",
         title = "Release by employment, citizenship and police checks",
         subtitle = "Dashed line is the overall release rate. Point size shows the number of arrests")
}

plot_colour_within_groups <- function() {
  cec <- rel_by(A, colour, employed, citizen) %>%
    mutate(group = paste0("Employed: ", employed, "\nCitizen: ", citizen),
           group = factor(group, levels = unique(group[order(p)])))
  ggplot(cec, aes(x = p, y = group, colour = colour)) +
    geom_line(aes(group = group), colour = "grey70", linewidth = 1.2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, alpha = 0.6) +
    geom_point(aes(size = n)) +
    scale_colour_manual(values = col_race, name = "Colour") +
    scale_size_area(max_size = 7, name = "Arrests") +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Percent released with a summons (95% CI)", y = NULL,
         title = "Release by colour within employment and citizenship groups")
}

plot_doubledecker <- function() {
  doubledecker(released ~ citizen + employed + colour, data = A,
               gp = gpar(fill = col_rel),
               labeling_args = list(abbreviate_labs = c(colour = 1)),
               main = "Release by citizenship, employment and colour (B = Black, W = White)")
}

plot_mosaic <- function() {
  mosaic(~ colour + employed + citizen + released, data = A,
         shade = TRUE, legend = TRUE,
         labeling_args = list(abbreviate_labs = c(colour = FALSE, employed = 3,
                                                  citizen = 3, released = 3)),
         main = "Mosaic of colour, employment, citizenship and release")
}

plot_checks_by_colour <- function() {
  A %>%
    count(colour, checks) %>%
    group_by(colour) %>%
    mutate(p = n / sum(n)) %>%
    ggplot(aes(x = factor(checks), y = p, fill = colour)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = col_race, name = "Colour") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Number of databases with the person's name", y = "Percent of arrestees in group",
         title = "Distribution of police checks by colour")
}

plot_forest <- function() {
  tidy_or <- function(fit, model) {
    est <- coef(fit)[-1]
    se  <- sqrt(diag(vcov(fit)))[-1]
    data.frame(term = names(est), OR = exp(est),
               lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se), model = model)
  }
  labs_map <- c(colourWhite = "White (vs Black)", employedYes = "Employed",
                citizenYes = "Citizen", checks = "Each additional check",
                age = "Age (per year)", sexMale = "Male", year = "Year (per year)")
  fo <- tidy_or(fit_adj, "Adjusted") %>%
    bind_rows(tidy_or(fit_crude, "Crude")) %>%
    mutate(term = factor(labs_map[term], levels = rev(labs_map)))
  ggplot(fo, aes(x = OR, y = term, colour = model)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2, 3)) +
    scale_colour_manual(values = c(Adjusted = "#2F4858", Crude = "#B0B0B0"), name = NULL) +
    labs(x = "Odds ratio for release with a summons (log scale, 95% CI)", y = NULL,
         title = "Logistic regression of release on all seven variables")
}

plot_age_period <- function() {
  rel_by(A, colour, agegrp, period) %>%
    ggplot(aes(x = agegrp, y = p, colour = colour, group = colour)) +
    geom_line(linewidth = 1) +
    geom_point(aes(size = n)) +
    facet_wrap(~ period) +
    scale_colour_manual(values = col_race, name = "Colour") +
    scale_size_area(max_size = 5, name = "Arrests") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Age group", y = "Percent released",
         title = "Release by colour across age groups and periods")
}

plot_distributions <- function() {
  p_fac <- A %>%
    select(colour, sex, employed, citizen) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "level") %>%
    count(variable, level) %>%
    group_by(variable) %>%
    mutate(p = n / sum(n)) %>%
    ungroup() %>%
    mutate(variable = factor(variable, levels = c("colour", "sex", "employed", "citizen"))) %>%
    ggplot(aes(x = level, y = n)) +
    geom_col(fill = "grey45", width = 0.65) +
    geom_text(aes(label = pct(p)), vjust = -0.3, size = 3.6) +
    facet_wrap(~ variable, nrow = 1, scales = "free_x") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Arrests", title = "Binary characteristics")
  p_age <- ggplot(A, aes(x = age)) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "grey45", colour = "white") +
    geom_vline(xintercept = median(A$age), linetype = 2) +
    labs(x = "Age (years)", y = "Arrests", title = "Age (1 year bins, dashed line = median)")
  p_chk <- ggplot(A, aes(x = factor(checks))) +
    geom_bar(fill = "grey45", width = 0.7) +
    labs(x = "Number of police checks", y = "Arrests", title = "Checks")
  p_yr <- ggplot(A, aes(x = factor(year))) +
    geom_bar(fill = "grey45", width = 0.7) +
    labs(x = "Year", y = "Arrests", title = "Year")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
  print(p_fac, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  print(p_age, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(p_chk, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p_yr,  vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2:3))
}

plot_spines <- function() {
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), cex.main = 1.1)
  on.exit(par(op))
  spineplot(released ~ colour,         data = A, col = col_rel, main = "Colour")
  spineplot(released ~ employed,       data = A, col = col_rel, main = "Employed")
  spineplot(released ~ citizen,        data = A, col = col_rel, main = "Citizen")
  spineplot(released ~ sex,            data = A, col = col_rel, main = "Sex")
  spineplot(released ~ factor(checks), data = A, col = col_rel, main = "Checks", xlab = "checks")
  spineplot(released ~ factor(year),   data = A, col = col_rel, main = "Year", xlab = "year")
}

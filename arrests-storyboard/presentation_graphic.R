# One presentation graphic: the colour gap in release by age and period.
# Observed release rates with 95% Wilson CIs. The colour x age and colour x year
# interactions are tested in the lowest-AIC logistic model (mod4).
library(carData)
library(ggplot2)
library(dplyr)
data(Arrests)

mod1 <- glm(released ~ colour + year + age + sex + employed + citizen + checks,
            family = binomial, data = Arrests)
mod4 <- glm(released ~ colour + year + age + employed + citizen + checks +
              colour:year + colour:age, family = binomial, data = Arrests)
lrt <- anova(update(mod4, . ~ . - colour:year - colour:age), mod4, test = "Chisq")

col_race <- c(Black = "#3B6FB6", White = "#E1A33B")

d <- Arrests %>%
  mutate(agegrp = cut(age, c(11, 17, 20, 25, 30, 40, 70),
                      labels = c("12-17", "18-20", "21-25", "26-30", "31-40", "41+")),
         period = factor(ifelse(year <= 1999, "1997 to 1999", "2000 to 2002"))) %>%
  group_by(period, colour, agegrp) %>%
  summarise(n = n(), k = sum(released == "Yes"), .groups = "drop") %>%
  mutate(p = k / n, z = 1.96,
         lo = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
         hi = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
         x = as.numeric(agegrp) + ifelse(colour == "Black", -0.08, 0.08))

# Gap annotation at ages 21-25 in each period
ann <- d %>%
  filter(agegrp == "21-25") %>%
  select(period, colour, p) %>%
  tidyr::pivot_wider(names_from = colour, values_from = p) %>%
  mutate(x = 3, label = sprintf("%.0f point\ngap", 100 * (White - Black)))

ends <- d %>% filter(agegrp == "41+")

g <- ggplot(d, aes(x = x, y = p, colour = colour)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.7, alpha = 0.55) +
  geom_line(linewidth = 1.6) +
  geom_point(size = 3.6) +
  geom_segment(data = ann, aes(x = x - 0.32, xend = x - 0.32, y = Black, yend = White),
               inherit.aes = FALSE, colour = "grey25", linewidth = 0.8,
               arrow = arrow(ends = "both", length = unit(0.13, "cm"))) +
  geom_text(data = ann, aes(x = x - 0.42, y = (Black + White) / 2, label = label),
            inherit.aes = FALSE, hjust = 1, size = 5.2, fontface = "bold",
            colour = "grey15", lineheight = 0.95) +
  geom_text(data = ends, aes(x = x + 0.18, label = colour), hjust = 0,
            size = 5.5, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ period) +
  scale_colour_manual(values = col_race, guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = levels(d$agegrp), limits = c(0.4, 7.1)) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0.5, 1, 0.1)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "Young Black arrestees were much less likely to be released with a summons",
       subtitle = "The gap closed by about age 40, and narrowed after 1999",
       x = "Age of arrestee (years)",
       y = "Percent released with a summons",
       caption = paste0(
         "Observed percent released with 95% Wilson CIs. Toronto marijuana possession arrests, 1997 to 2002, n = 5,226 (carData::Arrests).\n",
         "In the lowest-AIC logistic model (adjusting for year, age, employment, citizenship and police checks), colour x age and colour x year\n",
         sprintf("interactions improve fit (likelihood ratio chi-square = %.1f, 2 df, p < 0.0001). AIC %.0f vs %.0f without them.",
                 lrt$Deviance[2], AIC(mod4), AIC(update(mod4, . ~ . - colour:year - colour:age))))) +
  theme_minimal(base_size = 17) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(2.5, "lines"),
        strip.text = element_text(size = 17, face = "bold", hjust = 0),
        plot.title = element_text(face = "bold", size = 21),
        plot.subtitle = element_text(size = 16, colour = "grey30", margin = margin(b = 12)),
        plot.caption = element_text(size = 10.5, colour = "grey40", hjust = 0, lineheight = 1.1),
        plot.title.position = "plot", plot.caption.position = "plot")

ggsave("presentation_graphic.png", g, width = 13.33, height = 7.5, dpi = 150, bg = "white")

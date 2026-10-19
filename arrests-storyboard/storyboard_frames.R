# Export story frames 4 to 6 sized for a half-page storyboard box (7 x 4.6 in).
source("_setup.R")

box_theme <- theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(face = "bold", size = 17),
        legend.position = "top", legend.justification = "left")

save_box <- function(p, file, h = 4.6) ggsave(file, p, width = 7, height = h, dpi = 200, bg = "white")

# Frame 4: colour gap within employment x citizenship groups
cec <- rel_by(A, colour, employed, citizen) %>%
  mutate(group = paste0(ifelse(employed == "Yes", "Employed", "Unemployed"), ", ",
                        ifelse(citizen == "Yes", "citizen", "non-citizen")),
         group = factor(group, levels = unique(group[order(p)])))
gaps <- cec %>% select(group, colour, p) %>%
  tidyr::pivot_wider(names_from = colour, values_from = p) %>%
  mutate(lab = sprintf("%.0f pts", 100 * (White - Black)))
f4 <- ggplot(cec, aes(x = p, y = group)) +
  geom_line(aes(group = group), colour = "grey75", linewidth = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = colour), height = 0.15, alpha = 0.5) +
  geom_point(aes(colour = colour, size = n)) +
  geom_text(data = gaps, aes(x = (Black + White) / 2, y = group, label = lab),
            nudge_y = 0.32, size = 4.3, colour = "grey20") +
  scale_colour_manual(values = col_race, name = NULL) +
  scale_size_area(max_size = 7, guide = "none") +
  scale_x_continuous(labels = scales::percent, limits = c(0.4, 1)) +
  labs(x = "Percent released with a summons (95% CI)", y = NULL,
       title = "Release by colour within subgroups") +
  box_theme
save_box(f4, "frame4.png")

# Frame 5: crude vs adjusted odds ratios
f5 <- plot_forest() +
  labs(title = "Crude and adjusted odds ratios for release",
       x = "Odds ratio (log scale, 95% CI)") +
  box_theme
save_box(f5, "frame5.png")

# Frame 6: the presentation graphic, resized for the box
source("presentation_graphic.R")
f6 <- g + labs(title = "Release by colour, age and period", subtitle = NULL,
               caption = "Observed percent released with 95% CIs.") +
  theme(plot.title = element_text(face = "bold", size = 17),
        axis.title = element_text(size = 14), axis.text = element_text(size = 11.5),
        strip.text = element_text(size = 14, face = "bold", hjust = 0),
        plot.caption = element_text(size = 10, colour = "grey40", hjust = 0),
        panel.spacing = unit(1.2, "lines"),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
f6$layers[[5]]$aes_params$size <- 4    # gap labels
f6$layers[[6]]$aes_params$size <- 4.3  # Black / White end labels
save_box(f6, "frame6.png", h = 4.8)

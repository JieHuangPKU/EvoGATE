suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

metrics_path <- "results/FigureS4C_threshold_sweep/tables/threshold_sweep_metrics.tsv"
optima_path <- "results/FigureS4C_threshold_sweep/tables/threshold_optima.tsv"
output_path <- "results/FigureS4C_threshold_sweep/plots/FigureS4C_threshold_sweep.pdf"

metrics_df <- read_tsv(metrics_path, show_col_types = FALSE) %>%
  mutate(
    metric = factor(metric, levels = c("MCC", "F1")),
    fusion_method = factor(
      fusion_method,
      levels = c("Concat", "Gated", "Gated+WBCE", "Residual gated", "Residual gated+WBCE")
    )
  ) %>%
  filter(metric %in% c("MCC", "F1"))

optima_df <- read_tsv(optima_path, show_col_types = FALSE) %>%
  filter(threshold_type %in% c("val_f1_opt", "val_mcc_opt")) %>%
  mutate(
    marker_label = recode(
      threshold_type,
      val_f1_opt = "Validation F1 optimum",
      val_mcc_opt = "Validation MCC optimum"
    ),
    marker_label = factor(marker_label, levels = c("Validation F1 optimum", "Validation MCC optimum")),
    fusion_method = factor(
      fusion_method,
      levels = c("Concat", "Gated", "Gated+WBCE", "Residual gated", "Residual gated+WBCE")
    )
  )

marker_df <- metrics_df %>%
  inner_join(optima_df, by = c("fusion_method", "threshold"), relationship = "many-to-many")

palette_values <- c(
  "Concat" = "#3b5b92",
  "Gated" = "#2f7d4a",
  "Gated+WBCE" = "#c26d1a",
  "Residual gated" = "#8a4f9e",
  "Residual gated+WBCE" = "#b23a48"
)

shape_values <- c(
  "Validation F1 optimum" = 17,
  "Validation MCC optimum" = 15
)

plot_obj <- ggplot(metrics_df, aes(x = threshold, y = mean, color = fusion_method, fill = fusion_method)) +
  geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey35", linewidth = 0.5) +
  geom_point(
    data = marker_df,
    aes(shape = marker_label),
    size = 2.2,
    stroke = 0.5
  ) +
  facet_wrap(~metric, nrow = 1, scales = "fixed") +
  scale_color_manual(values = palette_values, drop = FALSE) +
  scale_fill_manual(values = palette_values, drop = FALSE) +
  scale_shape_manual(values = shape_values, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "Figure S4C. Threshold sensitivity of decision metrics",
    subtitle = "MCC and F1 vary across classification thresholds; AUROC/AUPRC are threshold-invariant ranking metrics.",
    x = "Decision threshold",
    y = "Performance",
    color = "Fusion method",
    fill = "Fusion method",
    shape = "Optimized threshold marker"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank()
  )

ggsave(output_path, plot_obj, width = 11, height = 4.8, units = "in", device = cairo_pdf)

cat(normalizePath(output_path), "\n")

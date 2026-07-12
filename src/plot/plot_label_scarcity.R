#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option("--summary-dir", type = "character"),
  make_option("--output-dir", type = "character")
)
args <- parse_args(OptionParser(option_list = option_list))

summary_dir <- normalizePath(args$`summary-dir`, mustWork = TRUE)
output_dir <- args$`output-dir`
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- read_tsv(file.path(summary_dir, "label_scarcity_summary.tsv"), show_col_types = FALSE)
ranking_df <- read_tsv(file.path(summary_dir, "label_scarcity_ranking_table.tsv"), show_col_types = FALSE)

model_levels <- c("GraphSAGE", "MLP", "SVM", "RF", "N2V_MLP", "DC", "CC")
palette_map <- c(
  "GraphSAGE" = "#1296BA",
  "MLP" = "#6F6F6F",
  "SVM" = "#8C8C8C",
  "RF" = "#4F4F4F",
  "N2V_MLP" = "#8EC7E8",
  "DC" = "#E8B36A",
  "CC" = "#B9B9B9"
)

plot_df <- summary_df %>%
  filter(model %in% model_levels) %>%
  mutate(
    model = factor(model, levels = model_levels),
    train_fraction_pct = train_fraction * 100
  )

line_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

save_plot <- function(plot_obj, file_stub, width = 7.2, height = 4.6) {
  ggsave(
    filename = file.path(output_dir, paste0(file_stub, ".pdf")),
    plot = plot_obj,
    device = cairo_pdf,
    width = width,
    height = height
  )
  ggsave(
    filename = file.path(output_dir, paste0(file_stub, ".png")),
    plot = plot_obj,
    dpi = 300,
    width = width,
    height = height
  )
}

auprc_plot <- ggplot(plot_df, aes(x = train_fraction_pct, y = mean_AUPRC, color = model, fill = model, group = model)) +
  geom_ribbon(aes(ymin = pmax(mean_AUPRC - sd_AUPRC, 0), ymax = pmin(mean_AUPRC + sd_AUPRC, 1)), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.0) +
  scale_color_manual(values = palette_map) +
  scale_fill_manual(values = palette_map) +
  scale_x_continuous(breaks = seq(10, 90, by = 10), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01), limits = c(0, NA)) +
  labs(
    title = "Label Scarcity Benchmark: AUPRC",
    x = "Training fraction",
    y = "AUPRC"
  ) +
  line_theme

auroc_plot <- ggplot(plot_df, aes(x = train_fraction_pct, y = mean_AUROC, color = model, fill = model, group = model)) +
  geom_ribbon(aes(ymin = pmax(mean_AUROC - sd_AUROC, 0), ymax = pmin(mean_AUROC + sd_AUROC, 1)), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.0) +
  scale_color_manual(values = palette_map) +
  scale_fill_manual(values = palette_map) +
  scale_x_continuous(breaks = seq(10, 90, by = 10), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
  labs(
    title = "Label Scarcity Benchmark: AUROC",
    x = "Training fraction",
    y = "AUROC"
  ) +
  line_theme

retention_plot <- ggplot(plot_df, aes(x = train_fraction_pct, y = performance_retention_AUPRC, color = model, fill = model, group = model)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#9A9A9A", linewidth = 0.5) +
  geom_ribbon(aes(ymin = pmax(performance_retention_AUPRC - (sd_AUPRC / baseline_AUPRC_at_0.90), 0), ymax = performance_retention_AUPRC + (sd_AUPRC / baseline_AUPRC_at_0.90)), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.0) +
  scale_color_manual(values = palette_map) +
  scale_fill_manual(values = palette_map) +
  scale_x_continuous(breaks = seq(10, 90, by = 10), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01), limits = c(0, NA)) +
  labs(
    title = "Label Scarcity Benchmark: AUPRC Retention",
    x = "Training fraction",
    y = "AUPRC retention vs 90%"
  ) +
  line_theme

save_plot(auprc_plot, "label_scarcity_auprc")
save_plot(auroc_plot, "label_scarcity_auroc")
save_plot(retention_plot, "label_scarcity_retention")

write_tsv(ranking_df, file.path(output_dir, "label_scarcity_ranking_table.tsv"))

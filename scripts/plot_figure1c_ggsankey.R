#!/usr/bin/env Rscript
# 脚本名称: plot_figure1c_ggsankey.R
# 日期: 2026-05-06
# 作者: OpenAI/Codex + Jie Huang workflow support
# 功能描述: 使用 ggsankey 绘制 Fusarium new label source-resolved transfer Sankey
#           Figure 1C，并在右侧添加 summary / scale box。

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsankey)
  library(patchwork)
  library(grid)
  library(svglite)
})

# 读取命令行参数并定位项目路径
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  script_path <- "scripts/plot_figure1c_ggsankey.R"
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
result_dir <- file.path(repo_root, "results", "Figure1C")

sankey_tsv <- file.path(result_dir, "figure1c_sankey_long.tsv")
summary_tsv <- file.path(result_dir, "figure1c_summary_box.tsv")
stage_tsv <- file.path(result_dir, "figure1c_stage_counts.tsv")

pdf_out <- file.path(result_dir, "Figure1C_ggsankey.pdf")
svg_out <- file.path(result_dir, "Figure1C_ggsankey.svg")
png_out <- file.path(result_dir, "Figure1C_ggsankey.png")

# 检查输入文件
required_files <- c(sankey_tsv, summary_tsv, stage_tsv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing Figure 1C plotting inputs:\n", paste(missing_files, collapse = "\n"))
}

# 读取 Python 整理后的 Sankey 与 summary 数据
sankey <- read.delim(sankey_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
summary_box <- read.delim(summary_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stage_counts <- read.delim(stage_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# 整理因子顺序与节点标签
stage_levels <- c("Support source", "Supported orthogroups", "Mapped Fusarium genes", "Final label class")
source_levels <- c(
  "S. cerevisiae only",
  "S. pombe only",
  "Shared by both yeasts",
  "PHI-base essential/lethal"
)
sankey$x <- factor(sankey$x, levels = stage_levels)
sankey$next_x <- factor(sankey$next_x, levels = stage_levels)
sankey$support_source <- factor(sankey$support_source, levels = source_levels)

# 构造稳定的柔和配色
source_palette <- c(
  "S. cerevisiae only" = "#5A8BBF",
  "S. pombe only" = "#6AA879",
  "Shared by both yeasts" = "#8B79B7",
  "PHI-base essential/lethal" = "#E2A84D"
)
final_palette <- c(
  "Essential (positive)" = "#A8DDB5",
  "Non-essential / excluded" = "#D6D6D6"
)

# 节点使用 source 色，终点节点使用 final class 色
sankey$node_fill <- as.character(sankey$support_source)
sankey$next_node_fill <- ifelse(
  sankey$next_x == "Final label class",
  sankey$next_node,
  as.character(sankey$support_source)
)
node_palette <- c(source_palette, final_palette)

# 准备层标题位置
layer_titles <- data.frame(
  x = factor(stage_levels, levels = stage_levels),
  label = stage_levels,
  y = c(1118, 1118, 1118, 1118)
)

# 绘制 Figure 1C Sankey 主体
sankey_plot <- ggplot(
  sankey,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    value = value,
    fill = support_source
  )
) +
  geom_sankey(
    width = 0.12,
    smooth = 7,
    alpha = 0.58,
    color = NA,
    space = 22
  ) +
  geom_sankey_label(
    aes(label = node, fill = node_fill),
    width = 0.12,
    space = 22,
    color = "#333333",
    size = 2.35,
    label.r = unit(0.02, "lines"),
    label.padding = unit(0.12, "lines"),
    family = "sans",
    show.legend = FALSE
  ) +
  geom_text(
    data = layer_titles,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold",
    color = "#222222"
  ) +
  scale_fill_manual(values = node_palette, drop = FALSE) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.16))) +
  labs(
    title = "Figure 1C. Source-resolved transfer path of Fusarium new positives",
    subtitle = "Positive labels are traced from support source through supported orthogroups and canonical Fusarium mapping to final label class.",
    fill = "Support source"
  ) +
  theme_sankey(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.4),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 13, face = "bold", color = "#1F1F1F", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.2, color = "#4B4B4B", margin = margin(b = 10)),
    plot.margin = margin(8, 12, 8, 24)
  )

# 整理 summary / scale box 文本
summary_order <- c(
  "old positives",
  "new positives",
  "high-confidence transferred positives",
  "old total",
  "new total",
  "PHI-supported positives",
  "yeast-transfer-supported positives"
)
summary_box$metric <- factor(summary_box$metric, levels = summary_order)
summary_box <- summary_box[order(summary_box$metric), ]
metric_labels <- c(
  "old positives" = "old positives",
  "new positives" = "new positives",
  "high-confidence transferred positives" = "HC transferred positives",
  "old total" = "old total",
  "new total" = "new total",
  "PHI-supported positives" = "PHI-supported positives",
  "yeast-transfer-supported positives" = "yeast-transfer positives"
)
summary_box$display_metric <- metric_labels[as.character(summary_box$metric)]
summary_box$label <- sprintf("%s  %s", summary_box$display_metric, format(summary_box$count, big.mark = ",", trim = TRUE))
summary_box$y <- rev(seq_len(nrow(summary_box)))

positive_sources <- subset(stage_counts, layer == "Support source")
positive_sources <- positive_sources[positive_sources$node %in% source_levels, ]
positive_sources$node <- factor(positive_sources$node, levels = source_levels)
positive_sources <- positive_sources[order(positive_sources$node), ]
positive_sources$source_label <- sprintf("%s: %s", positive_sources$node, positive_sources$count)
positive_sources$y <- rev(seq_len(nrow(positive_sources))) - 0.2

# 绘制右侧 summary / scale box
summary_plot <- ggplot() +
  annotate(
    "rect",
    xmin = 0.05,
    xmax = 0.95,
    ymin = 0.1,
    ymax = 9.35,
    fill = "#FBFBFB",
    color = "#8E8E8E",
    linewidth = 0.45,
    linetype = "22"
  ) +
  annotate(
    "text",
    x = 0.1,
    y = 8.95,
    label = "Summary / scale",
    hjust = 0,
    size = 3.6,
    fontface = "bold",
    color = "#222222"
  ) +
  geom_text(
    data = summary_box,
    aes(x = 0.1, y = y + 1.1, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.9,
    color = "#333333"
  ) +
  annotate(
    "segment",
    x = 0.1,
    xend = 0.9,
    y = 2.05,
    yend = 2.05,
    color = "#C8C8C8",
    linewidth = 0.35
  ) +
  annotate(
    "text",
    x = 0.1,
    y = 1.7,
    label = "Source composition",
    hjust = 0,
    size = 3.1,
    fontface = "bold",
    color = "#222222"
  ) +
  geom_text(
    data = positive_sources,
    aes(x = 0.1, y = y * 0.28 + 0.2, label = source_label, color = node),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.5
  ) +
  scale_color_manual(values = source_palette, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 9.5), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(44, 8, 34, 4)
  )

# 拼接主图和右侧 summary box
figure1c <- sankey_plot + summary_plot + plot_layout(widths = c(4.55, 1.45))

# 导出 PDF、SVG 和 PNG
ggsave(pdf_out, figure1c, width = 13.2, height = 7.1, device = cairo_pdf, bg = "white")
svglite(svg_out, width = 13.2, height = 7.1, bg = "white")
print(figure1c)
dev.off()
ggsave(png_out, figure1c, width = 13.2, height = 7.1, dpi = 360, bg = "white")

message("Figure 1C plot exported to ", pdf_out)
message("Figure 1C plot exported to ", svg_out)
message("Figure 1C plot exported to ", png_out)
message("Done.")

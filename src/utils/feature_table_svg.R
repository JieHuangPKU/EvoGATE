#!/usr/bin/env Rscript

# 加载必要的库
if (!require("readr")) install.packages("readr")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("grid")) install.packages("grid")
if (!require("svglite")) install.packages("svglite")

library(readr)
library(gridExtra)
library(grid)
library(svglite)

# 设置输出文件
output_dir <- "tables"
# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
svg_file1 <- file.path(output_dir, "feature_stats_table.svg")
svg_file2 <- file.path(output_dir, "organism_stats_table.svg")

# 检查文件是否存在
feature_file <- "outputs/results/feature_stats_20250420_092140.csv"
organism_file <- "outputs/results/organism_stats_20250420_085017.csv"

if (!file.exists(feature_file)) {
  stop(paste("文件不存在:", feature_file))
}
if (!file.exists(organism_file)) {
  stop(paste("文件不存在:", organism_file))
}

# 读取数据文件
tryCatch({
  feature_stats <- read.csv(feature_file, stringsAsFactors = FALSE)
  organism_stats <- read.csv(organism_file, stringsAsFactors = FALSE)
}, error = function(e) {
  stop(paste("读取文件错误:", e$message))
})

# 检查数据是否正确读取
if (nrow(feature_stats) == 0 || ncol(feature_stats) == 0) {
  stop("feature_stats 数据为空")
}
if (nrow(organism_stats) == 0 || ncol(organism_stats) == 0) {
  stop("organism_stats 数据为空")
}

# 重命名行名列
names(feature_stats)[1] <- "特征类型"
names(organism_stats)[1] <- "统计项"

# 格式化大数字
format_number <- function(df) {
  for(i in 2:ncol(df)) {
    if(is.numeric(df[,i])) {
      df[,i] <- format(df[,i], big.mark = ",", scientific = FALSE)
    }
  }
  return(df)
}

organism_stats <- format_number(organism_stats)

# 创建表格函数
create_table_grob <- function(df, title, colors) {
  # 表格背景色设置
  rows <- nrow(df)
  cols <- ncol(df)
  
  # 创建行和列的背景色矩阵
  bg_rows <- matrix(c(rep(c("white", "gray95"), length.out = rows)), nrow = rows, ncol = cols)
  bg_cols <- matrix(rep("white", rows * cols), nrow = rows, ncol = cols)
  
  # 设置表头背景色
  header_colors <- matrix(rep(colors$header, cols), nrow = 1, ncol = cols)
  
  # 设置第一列背景色
  first_col_colors <- matrix(rep(colors$first_col, rows), nrow = rows, ncol = 1)
  
  # 合并背景色
  bg_matrix <- rbind(header_colors, bg_rows)
  
  # 设置文本颜色
  text_colors <- matrix(rep("black", (rows+1) * cols), nrow = rows+1, ncol = cols)
  text_colors[1,] <- "white" # 表头文字为白色
  
  # 添加表头
  df_with_header <- rbind(names(df), df)
  
  # 创建表格grob
  grob <- gridExtra::tableGrob(
    df_with_header, 
    rows = NULL,
    theme = ttheme_minimal(
      core = list(
        bg_params = list(fill = bg_matrix, col = NA),
        fg_params = list(col = text_colors, fontface = c(rep("bold", cols), rep("plain", rows * cols)))
      ),
      colhead = list(fg_params = list(col = "white", fontface = "bold"))
    )
  )
  
  # 添加标题
  title_grob <- textGrob(title, gp = gpar(fontface = "bold", fontsize = 14))
  padding <- unit(0.5, "line")
  
  # 合并标题和表格
  result <- gtable::gtable_add_rows(
    grob, 
    heights = grobHeight(title_grob) + padding,
    pos = 0
  )
  result <- gtable::gtable_add_grob(
    result,
    title_grob,
    t = 1, l = 1, r = ncol(result)
  )
  
  return(result)
}

# 生成表格1
tryCatch({
  svglite(svg_file1, width = 10, height = 4)
  table1 <- create_table_grob(
    feature_stats, 
    "表1：各生物体特征维度统计", 
    list(header = "#4E79A7", first_col = "#A0CBE8")
  )
  grid.newpage()
  grid.draw(table1)
  dev.off()
  
  cat("特征维度表已保存为SVG:", svg_file1, "\n")
}, error = function(e) {
  cat("生成特征维度表失败:", e$message, "\n")
})

# 生成表格2
tryCatch({
  svglite(svg_file2, width = 10, height = 7)
  table2 <- create_table_grob(
    organism_stats, 
    "表2：各生物体网络和基因统计", 
    list(header = "#E15759", first_col = "#F1A8A7")
  )
  grid.newpage()
  grid.draw(table2)
  dev.off()
  
  cat("生物体统计表已保存为SVG:", svg_file2, "\n")
}, error = function(e) {
  cat("生成生物体统计表失败:", e$message, "\n")
}) 
#!/usr/bin/env Rscript

# 加载必要的库
if (!require("kableExtra")) install.packages("kableExtra")
if (!require("readr")) install.packages("readr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("knitr")) install.packages("knitr")
if (!require("viridis")) install.packages("viridis")

library(kableExtra)
library(readr)
library(dplyr)
library(knitr)
library(viridis)

# 设置工作目录，确保文件路径正确
# setwd("工作目录路径") # 如果需要，取消注释并设置适当的路径

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
  feature_stats <- read_csv(feature_file)
  organism_stats <- read_csv(organism_file)
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

# 设置输出目录
output_dir <- "tables"
# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 设置输出文件路径
html_file_path <- file.path(output_dir, "feature_tables.html")

# 创建HTML文件
html_file <- file(html_file_path, "w", encoding = "UTF-8")

# 写入HTML头部
cat('<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>生物体特征和网络统计表</title>
  <style>
    body { font-family: Arial, "Microsoft YaHei", sans-serif; margin: 20px; }
    .table-container { margin-bottom: 30px; }
    caption { font-weight: bold; font-size: 1.2em; margin-bottom: 10px; }
  </style>
</head>
<body>
<h1>生物体特征和网络统计表</h1>
<div class="table-container">
', file = html_file)

# 表1: 特征统计表
feature_table <- kbl(feature_stats, 
                    caption = "表1：各生物体特征维度统计",
                    format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#4E79A7") %>%
  add_header_above(c(" " = 1, "多组学特征维度" = ncol(feature_stats) - 1))

# 写入特征表格到HTML
cat(as.character(feature_table), file = html_file)
cat('</div>\n<div class="table-container">\n', file = html_file)

# 表2: 生物体统计表
# 格式化大数字
organism_stats[, 2:ncol(organism_stats)] <- lapply(organism_stats[, 2:ncol(organism_stats)], function(x) {
  format(x, big.mark = ",", scientific = FALSE)
})

organism_table <- kbl(organism_stats, 
                     caption = "表2：各生物体网络和基因统计",
                     format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#E15759") %>%
  add_header_above(c(" " = 1, "生物体网络统计" = ncol(organism_stats) - 1))

# 写入生物体表格到HTML
cat(as.character(organism_table), file = html_file)

# 写入HTML尾部
cat('
</div>
</body>
</html>', file = html_file)

# 关闭文件
close(html_file)

# 检查HTML是否成功生成
if (file.exists(html_file_path) && file.info(html_file_path)$size > 0) {
  cat("HTML文件已成功生成：", html_file_path, "\n")
  cat("请在浏览器中打开此文件查看表格，然后可以使用浏览器的打印功能导出为PDF\n")
} else {
  cat("HTML文件生成失败或文件为空\n")
} 
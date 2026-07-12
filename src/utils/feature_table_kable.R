#!/usr/bin/env Rscript

# 加载必要的库
if (!require("kableExtra")) install.packages("kableExtra")
if (!require("readr")) install.packages("readr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("knitr")) install.packages("knitr")
if (!require("viridis")) install.packages("viridis")
if (!require("rmarkdown")) install.packages("rmarkdown")

library(kableExtra)
library(readr)
library(dplyr)
library(knitr)
library(viridis)
library(rmarkdown)

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
temp_feature_file <- file.path(output_dir, "temp_feature_stats.csv")
temp_organism_file <- file.path(output_dir, "temp_organism_stats.csv")
rmd_file <- file.path(output_dir, "feature_tables_temp.Rmd")
output_file <- file.path(output_dir, "feature_tables.pdf")
html_output <- file.path(output_dir, "feature_tables.html")

write.csv(feature_stats, temp_feature_file, row.names = FALSE)
write.csv(organism_stats, temp_organism_file, row.names = FALSE)

# 创建临时Rmd文件
cat('---
title: "生物体特征和网络统计表"
output: 
  pdf_document:
    latex_engine: xelatex
    toc: false
    keep_tex: false
mainfont: Arial
CJKmainfont: Arial Unicode MS
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(kableExtra)
```

## 表1: 各生物体特征维度统计

```{r feature-table, echo=FALSE}
feature_stats <- read.csv("', temp_feature_file, '", stringsAsFactors = FALSE)

kbl(feature_stats, 
    booktabs = TRUE,
    align = rep("c", ncol(feature_stats))) %>%
  kable_styling(latex_options = c("striped"), 
                full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE) %>%
  add_header_above(c(" " = 1, "多组学特征维度" = ncol(feature_stats) - 1))
```

## 表2: 各生物体网络和基因统计

```{r organism-table, echo=FALSE}
organism_stats <- read.csv("', temp_organism_file, '", stringsAsFactors = FALSE)

# 格式化大数字
for(i in 2:ncol(organism_stats)) {
  if(is.numeric(organism_stats[,i])) {
    organism_stats[,i] <- format(organism_stats[,i], big.mark = ",", scientific = FALSE)
  }
}

kbl(organism_stats, 
    booktabs = TRUE,
    align = rep("c", ncol(organism_stats))) %>%
  kable_styling(latex_options = c("striped"), 
                full_width = FALSE) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE) %>%
  add_header_above(c(" " = 1, "生物体网络统计" = ncol(organism_stats) - 1))
```
', file = rmd_file)

# 使用rmarkdown渲染PDF
tryCatch({
  render(rmd_file, output_file = output_file)
  
  # 检查PDF是否成功生成
  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    cat("PDF已成功生成：", output_file, "\n")
  } else {
    cat("PDF生成失败或文件为空\n")
  }
  
  # 删除临时文件
  if (file.exists(rmd_file)) file.remove(rmd_file)
  if (file.exists(temp_feature_file)) file.remove(temp_feature_file)
  if (file.exists(temp_organism_file)) file.remove(temp_organism_file)
  
}, error = function(e) {
  cat("PDF生成错误:", e$message, "\n")
  cat("尝试使用HTML格式作为替代方案...\n")
  
  # 如果PDF失败，尝试生成HTML
  # 直接使用kableExtra生成HTML文件，不依赖于rmarkdown
  cat("创建HTML格式的表格...\n")
  
  # 创建HTML文件
  html_file <- file(html_output, "w", encoding = "UTF-8")
  
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
  organism_stats_formatted <- organism_stats
  for(i in 2:ncol(organism_stats_formatted)) {
    if(is.numeric(organism_stats_formatted[,i])) {
      organism_stats_formatted[,i] <- format(organism_stats_formatted[,i], big.mark = ",", scientific = FALSE)
    }
  }
  
  organism_table <- kbl(organism_stats_formatted, 
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
  
  # 提示信息
  if (file.exists(html_output) && file.info(html_output)$size > 0) {
    cat("HTML文件已成功生成：", html_output, "\n")
    cat("请在浏览器中打开此文件查看表格，然后可以使用浏览器的打印功能导出为PDF\n")
  } else {
    cat("HTML文件生成失败\n")
  }
  
  # 删除临时文件
  if (file.exists(rmd_file)) file.remove(rmd_file)
  if (file.exists(temp_feature_file)) file.remove(temp_feature_file)
  if (file.exists(temp_organism_file)) file.remove(temp_organism_file)
}) 
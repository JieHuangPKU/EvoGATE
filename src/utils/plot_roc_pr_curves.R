#!/usr/bin/env Rscript

# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(gridExtra)
library(cowplot)  # 显式导入cowplot包

# 设置工作目录
base_dir <- "/Users/jiehuang/work/2025禾谷镰刀菌/程序/ProGATE/results/2025"

# 按照training_epoches.R中的顺序排列物种列表
species_list <- c("fgraminearum", "scerevisiae", "ecoli", "human", "melanogaster", "celegans")

# 物种代码到学名的映射
species_names <- c(
  "fgraminearum" = "F. graminearum", 
  "scerevisiae" = "S. cerevisiae", 
  "ecoli" = "E. coli", 
  "human" = "H. sapiens", 
  "melanogaster" = "D. melanogaster", 
  "celegans" = "C. elegans"
)

# 需要处理的模型列表
model_list <- c("transformer", "lstm", "cnn", "esm", "gat")

# 使用指定的颜色映射
model_colors <- c(
  "transformer" = "#88CCEEFF",  # 浅蓝色
  "lstm" = "#CC6677FF",         # 红色
  "cnn" = "#999933FF",          # 黄绿色
  "esm" = "#117733FF",          # 绿色
  "gat" = "#332288FF"           # 深蓝色
)

# 创建一个空的数据框来存储所有汇总数据
all_data <- data.frame()

# 为每个物种和模型对预先设置默认的AUC和AUPR值
performance_data <- expand.grid(
  Species = species_list,
  Model = model_list,
  Mean_AUC = NA_real_,
  Mean_AUPR = NA_real_,
  stringsAsFactors = FALSE
)

# 首先确保Species和Model列是字符串类型
performance_data$Species <- as.character(performance_data$Species)
performance_data$Model <- as.character(performance_data$Model)

# 临时存储F. graminearum的数据，用于后续调整
fgram_data_list <- list()

# 遍历所有物种和模型
for (species in species_list) {
  for (model in model_list) {
    # 构建平均指标文件路径
    metrics_path <- file.path(base_dir, species, "results", model, "kfold_model_saving", "average_metrics.csv")
    performance_file <- file.path(base_dir, species, "results", model, "kfold_model_saving", "average_performance.txt")
    
    # 检查文件是否存在
    if (file.exists(metrics_path)) {
      tryCatch({
        # 读取CSV文件
        data <- read_csv(metrics_path)
        
        # 添加物种和模型信息列
        data$Species <- species
        data$Model <- model
        
        # 读取性能数据文件
        if (file.exists(performance_file)) {
          # 尝试读取性能数据
          perf_lines <- readLines(performance_file)
          auc_line <- grep("Average AUC", perf_lines, value = TRUE)
          aupr_line <- grep("Average AUPR", perf_lines, value = TRUE)
          
          if (length(auc_line) > 0 && length(aupr_line) > 0) {
            # 提取AUC和AUPR值
            auc_val <- as.numeric(gsub(".*: ([0-9.]+).*", "\\1", auc_line))
            aupr_val <- as.numeric(gsub(".*: ([0-9.]+).*", "\\1", aupr_line))
            
            # 将值添加到数据中
            data$Mean_AUC <- auc_val
            data$Mean_AUPR <- aupr_val
            
            # 更新性能数据
            idx <- which(performance_data$Species == species & performance_data$Model == model)
            if (length(idx) > 0) {
              performance_data$Mean_AUC[idx] <- auc_val
              performance_data$Mean_AUPR[idx] <- aupr_val
            }
          }
        }
        
        # 如果是F. graminearum，先临时存储
        if (species == "fgraminearum") {
          fgram_data_list[[model]] <- data
        } else {
          # 将数据添加到汇总数据框
          all_data <- bind_rows(all_data, data)
        }
        
        cat(sprintf("成功处理 %s 物种的 %s 模型数据\n", species, model))
      }, error = function(e) {
        cat(sprintf("处理 %s 物种的 %s 模型数据时发生错误: %s\n", species, model, e$message))
      })
    } else {
      cat(sprintf("找不到 %s 物种的 %s 模型平均指标文件\n", species, model))
    }
  }
}

# 调整F. graminearum的数据，将Mean_FPR乘以1.03，大于1的设为1
for (model in names(fgram_data_list)) {
  data <- fgram_data_list[[model]]
  if (!is.null(data)) {
    # 不做任何调整，直接添加到汇总数据框
    all_data <- bind_rows(all_data, data)
  }
}

# 绘制并保存图表
if (nrow(all_data) > 0) {
  # 确保字符列是字符类型
  all_data$Species <- as.character(all_data$Species)
  all_data$Model <- as.character(all_data$Model)
  
  # 如果没有Mean_Recall列，从Mean_TPR计算
  if (!"Mean_Recall" %in% colnames(all_data)) {
    all_data$Mean_Recall <- all_data$Mean_TPR
  }
  
  # 添加学名列
  all_data$Species_Name <- species_names[all_data$Species]
  all_data$Species_Name <- factor(all_data$Species_Name, 
                                levels = c("F. graminearum", "S. cerevisiae", "E. coli", 
                                          "H. sapiens", "D. melanogaster", "C. elegans"))
  
  # 绘制ROC曲线
  roc_plot <- ggplot(all_data, aes(x = Mean_FPR, y = Mean_TPR, color = Model)) +
    geom_line(linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    facet_wrap(~ Species_Name, ncol = 3) +
    scale_color_manual(values = model_colors) +
    labs(
      title = "ROC Curves for Different Species",
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14, color = "black"),
      legend.title = element_text(size = 10),
      strip.text = element_text(face = "bold.italic", size = 12),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 1)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  # 保存ROC曲线图
  ggsave("roc_curves_all_species.pdf", roc_plot, width = 11, height = 6, dpi = 300)
  cat("已生成所有物种的ROC曲线图\n")
  
  # 绘制PR曲线
  pr_plot <- ggplot(all_data, aes(x = Mean_Recall, y = Mean_Precision, color = Model)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ Species_Name, ncol = 3) +
    scale_color_manual(values = model_colors) +
    labs(
      title = "Precision-Recall Curves for Different Species",
      x = "Recall",
      y = "Precision"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14, color = "black"),
      legend.title = element_text(size = 10),
      strip.text = element_text(face = "bold.italic", size = 12),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 1)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  # 保存PR曲线图
  ggsave("pr_curves_all_species.pdf", pr_plot, width = 11, height = 6, dpi = 300)
  cat("已生成所有物种的PR曲线图\n")
  
  # 将汇总数据保存为CSV文件，以便进一步分析
  write_csv(all_data, "model_metrics_summary.csv")
  write_csv(performance_data, "model_performance_summary.csv")
  cat("已将汇总数据保存到 model_metrics_summary.csv 和 model_performance_summary.csv\n")
} else {
  cat("没有找到任何物种的模型数据，请检查文件路径和格式\n")
} 
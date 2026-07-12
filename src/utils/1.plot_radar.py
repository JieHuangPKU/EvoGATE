#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
雷达图生成脚本 - 从三个CSV文件(GAT_results_v2.csv, GCN_results_v2.csv, GIN_results_v2.csv)
绘制三个物种(yeast, human, celegans)下不同模型(GAT, GIN, GCN)的性能对比
每个子图比较4种特征：各模型的ESM2嵌入和组学特征(EXP_SUB_ORT)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.font_manager import FontProperties
from math import pi

# 设置PDF输出为可编辑模式
plt.rcParams["pdf.fonttype"] = 42  # 使用TrueType字体
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans", "sans-serif"]

# 创建Arial Bold斜体字体
arial_bold_italic = FontProperties(family="Arial", weight="bold", style="italic")

# 文件路径
GAT_CSV = "GAT_results_v2.csv"
GCN_CSV = "GCN_results_v2.csv"
GIN_CSV = "GIN_results_v2.csv"
OUTPUT_DIR = "plots"

# 确保输出目录存在
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 性能指标和物种名称
METRICS = ["AUROC", "AUCPR", "Precision", "MCC", "Accuracy"]
ORGANISMS = ["celegans", "yeast", "human"]
# 物种拉丁学名映射
ORGANISM_DISPLAY_NAMES = {
    "celegans": "C. elegans",
    "yeast": "S. cerevisiae",
    "human": "H. sapiens"
}
MODELS = ["GAT", "GIN", "GCN"]

# 颜色设置 - 专业配色方案
COLOR_MAP = {
    "GAT_ESM2": "#FF420EFF",      # 深红色
    "GIN_ESM2": "#579D1CFF",      # 深蓝色
    "GCN_ESM2": "#83CAFFFF",      # 紫色
    "GAT_EXP_SUB_ORT": "#4B1F6FFF", # 橙色
    "GIN_EXP_SUB_ORT": "#4B1F6FFF", # 绿色
    "GCN_EXP_SUB_ORT": "#4B1F6FFF"  # 青色

}

def load_data(csv_files):
    """加载并处理多个CSV数据"""
    try:
        dataframes = []
        for file in csv_files:
            df = pd.read_csv(file)
            dataframes.append(df)
        
        # 合并所有数据
        combined_df = pd.concat(dataframes, ignore_index=True)
        return combined_df
    except Exception as e:
        print(f"加载CSV文件时出错: {e}")
        return None

def create_radar_chart(ax, df, organism, metrics):
    """创建单个物种的雷达图"""
    # 选择该物种的数据
    org_data = df[df['organism'] == organism]
    
    # 特征列表
    features = []
    
    # 添加ESM2嵌入特征
    for model in MODELS:
        # 选择对应的ESM2嵌入模型行
        emb_data = org_data[org_data['name'] == f"{model}_ESM2"]
        if not emb_data.empty:
            features.append({'name': f"{model}_ESM2", 'data': emb_data})
    
    # 只添加GAT模型的组学特征数据
    omics_data = org_data[org_data['name'] == "GAT_EXP_SUB_ORT"]
    if not omics_data.empty:
        features.append({'name': "GAT_EXP_SUB_ORT", 'data': omics_data})
    
    # 计算角度
    N = len(metrics)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]  # 闭合图形
    
    # 设置第一个轴在顶部
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    
    # 绘制轴并添加标签
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics, size=14)
    
    # 设置y轴标签和限制
    ax.set_rlabel_position(0)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], color="grey", size=9)
    ax.set_ylim(0, 1.0)
    
    # 绘制每个特征的性能
    legend_handles = []
    
    for feature in features:
        feature_name = feature['name']
        feature_data = feature['data']
        
        if feature_data.empty:
            print(f"警告: 找不到{organism}的{feature_name}数据")
            continue
        
        # 获取性能指标
        values = []
        for metric in metrics:
            if metric == "AUROC":
                values.append(feature_data['mean'].iloc[0])
            elif metric == "AUCPR":
                values.append(feature_data['auc_pr'].iloc[0])
            elif metric == "Precision":
                values.append(feature_data['precision'].iloc[0])
            elif metric == "MCC":
                values.append(feature_data['mcc'].iloc[0])
            elif metric == "Accuracy":
                values.append(feature_data['accuracy'].iloc[0])
        
        # 闭合多边形
        values += values[:1]
        
        # 获取颜色
        color = COLOR_MAP.get(feature_name, "gray")
        
        # 绘制线和填充区域
        ax.plot(angles, values, linewidth=2, linestyle='solid', color=color, label=feature_name)
        ax.fill(angles, values, color=color, alpha=0.1)
        
        # 创建图例项
        legend_handles.append(Line2D([0], [0], color=color, lw=2, label=feature_name))
    
    # 添加标题 - 使用Arial Bold斜体字体
    display_name = ORGANISM_DISPLAY_NAMES.get(organism, organism.capitalize())
    ax.set_title(display_name, fontsize=12, fontproperties=arial_bold_italic, y=1.1)
    
    return legend_handles

def main():
    """主函数"""
    print(f"正在加载CSV数据文件...")
    df = load_data([GAT_CSV, GCN_CSV, GIN_CSV])
    
    if df is None:
        print("数据加载失败，退出程序")
        return
    
    print("正在创建雷达图...")
    
    # 创建1x3的图形布局
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), subplot_kw={'polar': True})
    
    # 所有图例句柄
    all_handles = []
    
    # 为每个物种创建雷达图
    for i, organism in enumerate(ORGANISMS):
        handles = create_radar_chart(axes[i], df, organism, METRICS)
        if i == 0:  # 只保存第一个物种的图例句柄
            all_handles = handles
    
    # 在图形底部添加共享图例
    fig.legend(handles=all_handles, loc='lower center', ncol=len(all_handles), 
               bbox_to_anchor=(0.5, 0), fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # 为底部的图例腾出空间
    
    # 添加总标题
    plt.suptitle('', fontsize=16, fontweight='bold', y=0.98)
    
    # 保存图形
    output_path = os.path.join(OUTPUT_DIR, "gnn_models_comparison_radar.pdf")
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    print(f"雷达图已保存至: {output_path}")
    
    plt.close()

if __name__ == "__main__":
    main()

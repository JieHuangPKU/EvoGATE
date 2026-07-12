import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

# 设置字体和风格
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# 设置所有文本为黑色
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.titlecolor'] = 'black'

# 确保输出目录存在
os.makedirs('outputs/plots', exist_ok=True)

# 读取数据，从各个物种的 _final.csv 文件
fg_df = pd.read_csv('outputs/results/fgraminearum_final.csv')
yeast_df = pd.read_csv('outputs/results/yeast_final.csv')
coli_df = pd.read_csv('outputs/results/coli_final.csv')
human_df = pd.read_csv('outputs/results/human_final.csv')
fly_df = pd.read_csv('outputs/results/melanogaster_final.csv')
celegans_df = pd.read_csv('outputs/results/celegans_final.csv')

# 定义要比较的方法
methods = {
    'GAT': 'GAT',
    'GCN': 'GCN',
    'DC': 'DC',      # 对应NDC开头
    'CC': 'CC'       # 对应NCC开头
}

# 定义要绘制的指标
metrics = {
    'ROC': 'mean',             # 使用mean表示ROC曲线下面积
    'AUPRC': 'auc_pr',         # PR曲线下面积
    'Specificity': 'precision', # 特异性
    'Accuracy': 'accuracy',    # 准确率
    'MCC': 'mcc'               # Matthews相关系数
}

# 设定标准差对应的列
std_metrics = {
    'ROC': 'std',
    'AUPRC': 'auc_pr_std',
    'Specificity': 'precision_std',
    'Accuracy': 'accuracy_std',
    'MCC': 'mcc_std'
}

# 设置参考图片中的颜色方案
# 根据参考图调整，红棕色(GAT)，深蓝色(GCN)，浅蓝色(DC)，浅棕色(CC)
bar_colors = ['#C15B4F', '#3B6C6D', '#A2A290', '#D8AE9C']

# 使用ggplot风格，但做一些自定义调整
plt.style.use('ggplot')
fig = plt.figure(figsize=(17, 10))

# 主图区域设置
ax = fig.add_subplot(111)
ax.set_facecolor('white')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax.set_ylabel('Performance', labelpad=15)

def plot(ax, df, title=None):
    """绘制一个子图，横坐标为不同指标，每个指标有不同模型的柱状图
    
    Args:
        ax: matplotlib Axes对象
        df: 包含结果的DataFrame
        title: 子图标题
    """
    # 设置柱状图的宽度和间距
    bar_width = 0.2
    x = np.arange(len(metrics))  # 指标位置
    
    # 设置背景样式 - 匹配参考图片
    ax.set_facecolor('#E5E5E5')  # 浅灰色背景
    ax.grid(color='white', linestyle='-', linewidth=1)  # 白色网格线
    ax.set_axisbelow(True)  # 确保网格线在柱状图下方
    
    # 保存所有方法的结果，用于绘图
    all_means = {}
    all_stds = {}
    
    # 对于每个方法，提取每个指标的均值和标准差
    for i, (method_key, method_name) in enumerate(methods.items()):
        means = []
        stds = []
        
        # 根据不同的方法类型处理数据获取
        if method_key == 'DC':
            # 查找NDC开头的行
            mask = df['name'].str.startswith('NDC_', na=False)
            method_df = df[mask]
            if not method_df.empty:
                method_df = method_df.iloc[[0]]  # 只取第一行
        elif method_key == 'CC':
            # 查找NCC开头的行
            mask = df['name'].str.startswith('NCC_', na=False)
            method_df = df[mask]
            if not method_df.empty:
                method_df = method_df.iloc[[0]]  # 只取第一行
        else:
            # 查找匹配的行（按方法名称）
            mask = (df['method'] == method_key)
            method_df = df[mask]
        
        if not method_df.empty:
            # 对于每个指标，提取均值和标准差
            for metric_name, column_name in metrics.items():
                std_column = std_metrics[metric_name]
                
                if column_name in method_df.columns and std_column in method_df.columns:
                    # 使用第一行数据
                    mean_val = method_df[column_name].values[0]
                    std_val = method_df[std_column].values[0]
                    
                    means.append(mean_val)
                    stds.append(std_val)
                else:
                    print(f"警告: 在{title}中没有找到方法{method_key}的{column_name}或{std_column}列")
                    means.append(np.nan)
                    stds.append(np.nan)
        else:
            print(f"警告: 在{title}中没有找到方法{method_key}的数据")
            means.extend([np.nan] * len(metrics))
            stds.extend([np.nan] * len(metrics))
        
        all_means[method_key] = means
        all_stds[method_key] = stds
    
    # 绘制柱状图
    for i, (method_key, method_name) in enumerate(methods.items()):
        # 计算柱状图的位置
        bar_positions = x - bar_width * 1.5 + bar_width * i
        
        # 绘制柱状图
        rects = ax.bar(bar_positions, all_means[method_key], width=bar_width, 
                     yerr=all_stds[method_key], align='center', 
                     color=bar_colors[i], capsize=3, label=method_key)
    
    # 设置坐标轴和标签
    ax.set_xticks(x)
    ax.set_xticklabels(list(metrics.keys()))
    ax.set_ylim(bottom=0.0, top=1.0)  # 设置y轴范围为0-1
    
    # 设置标题
    if title:
        ax.set_title(title, fontsize=13)

# 按照指定顺序排列物种
organisms = [
    ("F. graminearum", fg_df),
    ("S. cerevisiae", yeast_df),
    ("E. coli", coli_df),
    ("H. sapiens", human_df),
    ("D. melanogaster", fly_df),
    ("C. elegans", celegans_df)
]

# 创建 3x2 网格的子图
for i, (organism_name, df) in enumerate(organisms):
    ax = fig.add_subplot(2, 3, i+1)  # 2行3列
    plot(ax, df, organism_name)

# 添加图例 - 放在右上角，与参考图一致
handles, labels = ax.get_legend_handles_labels()
if handles:  # 确保有图例元素
    # 创建图例框，添加标题"Models"
    legend = fig.legend(handles, labels, title="Models", loc='upper right', 
              bbox_to_anchor=(0.95, 0.95))
    # 修改图例框背景为白色
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('lightgray')

# 使用 tight_layout 调整子图之间的间距
plt.tight_layout()
plt.suptitle('Performance Comparison of Network Methods across Species', fontsize=15)
plt.subplots_adjust(top=0.92, right=0.89)  # 为标题和图例腾出空间

# 保存图表
plt.savefig('outputs/plots/network_comparison.pdf')
plt.savefig('outputs/plots/network_comparison.png', dpi=300)
plt.close()

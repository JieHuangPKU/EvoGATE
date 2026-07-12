import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import os 
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap.umap_ as umap

# 添加项目根目录到导入路径
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, "../.."))
sys.path.append(project_root)

# 添加config目录到导入路径
config_dir = os.path.join(project_root, "config")
if config_dir not in sys.path:
    sys.path.insert(0, config_dir)

# 导入配置管理
from config.config_manager import get_species_config

from utils import * 
from data_loader import SeqDataset, collate_fn
from torch.utils.data import DataLoader 
from model import fc_layer, FGM

# 设置随机种子
torch.cuda.manual_seed(1029)
torch.manual_seed(1223)

def get_layer_output(model, input_data, layer_name):
    """获取模型中指定层的输出"""
    outputs = {}
    def get_output(layer, input, output):
        outputs[layer_name] = output.detach()  # 使用detach()避免构建计算图
    target_layer = dict([*model.named_modules()])[layer_name]
    hook = target_layer.register_forward_hook(get_output)
    model(input_data)
    hook.remove()
    return outputs[layer_name]

def transfer_list_tensor(b_x):
    """将列表张量转换为批量张量"""
    batch_tensor = torch.zeros((len(b_x), b_x[0].size(0)))  # [batch_size,1280]
    for row_id in range(len(b_x)):
        batch_tensor[row_id] = b_x[row_id]
    return batch_tensor

def extract_features(model, data_loader, device):
    """从数据中提取特征"""
    model.eval()
    layer_outputs = []
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    
    with torch.no_grad():
        for idx, (batch_data, batch_label) in enumerate(data_loader):
            batch_data = transfer_list_tensor(batch_data).to(device)
            batch_label = batch_label.float().to(device)
            output = model(batch_data)
            fc_g1_output = get_layer_output(model, batch_data, "fc_g1")
            layer_outputs.append(fc_g1_output.cpu().detach().numpy())
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat((total_labels, batch_label.view(-1, 1).cpu()), 0)
        
        total_labels_arr = total_labels.numpy().flatten()
        total_preds_arr = total_preds.detach().numpy().flatten() 

    features_array = np.vstack(layer_outputs)
    print(f"提取的特征形状: {features_array.shape}")
    
    return features_array, total_labels_arr

def visualize_and_save(features, labels, output_prefix, species):
    """生成并保存可视化图表"""
    # 使物种名称首字母大写，便于在标题中展示
    species_display = species.capitalize()
    
    # 限制样本数量以加快处理
    max_samples = min(10000, len(features))
    
    # 1. 直接t-SNE可视化
    print("执行t-SNE降维...")
    tsne = TSNE(n_components=2, random_state=42)
    tsne_results = tsne.fit_transform(features[:max_samples])

    plt.figure(figsize=(7, 7))
    plt.scatter(tsne_results[labels[:max_samples] == 0, 0], 
                tsne_results[labels[:max_samples] == 0, 1],
                s=1, c='tab:cyan', label='Non-essential')
    plt.scatter(tsne_results[labels[:max_samples] == 1, 0], 
                tsne_results[labels[:max_samples] == 1, 1],
                s=5, c='tab:red', label='Essential')
    plt.xlabel("t-SNE Dimension 1")
    plt.ylabel("t-SNE Dimension 2")
    plt.title(f"ESM model of {species_display} E/NE Genes")
    plt.legend()
    plt.savefig(f"{output_prefix}_t-SNE.pdf", format='pdf')
    plt.close()  # 关闭图形以释放内存

    # 2. UMAP可视化
    print("执行UMAP降维...")
    umap_model = umap.UMAP(n_components=2, verbose=True)
    umap_results = umap_model.fit_transform(features[:max_samples])

    plt.figure(figsize=(7, 7))
    plt.scatter(umap_results[labels[:max_samples] == 0, 0], 
                umap_results[labels[:max_samples] == 0, 1], 
                s=1, c='xkcd:azure', label='Non-essential')
    plt.scatter(umap_results[labels[:max_samples] == 1, 0], 
                umap_results[labels[:max_samples] == 1, 1], 
                s=5, c='xkcd:strong pink', label='Essential')
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")
    plt.title(f"ESM model of {species_display} E/NE Genes")
    plt.legend()
    plt.savefig(f"{output_prefix}_UMAP.pdf", format='pdf')
    plt.close()  # 关闭图形以释放内存
    
    print(f"可视化图表已保存到 {output_prefix}_*.pdf 文件")

def main():
    """主函数：使用argparse解析命令行参数"""
    # 使用argparse解析命令行参数
    parser = argparse.ArgumentParser(description='ESM特征提取与可视化')
    parser.add_argument('-species', '--species', type=str, required=True, help='物种名称')
    parser.add_argument('-fold', '--fold', type=int, default=0, help='指定使用的fold索引(默认为0)')
    parser.add_argument('-output_dir', '--output_dir', type=str, help='输出目录路径')
    parser.add_argument('-num_threads', '--num_threads', type=int, default=16, help='PyTorch内部线程数')
    
    # 使用parse_known_args()允许未知参数
    args, unknown = parser.parse_known_args()
    if unknown:
        print(f"警告: 忽略未知参数: {unknown}")
    
    # 设置设备
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type == 'cpu':
        torch.set_num_threads(args.num_threads)
        print(f"设置PyTorch线程数为: {args.num_threads}")
    
    # 获取配置
    config = get_species_config(args.species)
    print(f"成功加载 {args.species} 的配置")
    
    # 设置参数
    species = config.species_topofallfeature
    fold = args.fold
    kfold_root_path = config.kfold_root_path_topofallfeature
    
    # 设置模型保存路径
    model_saving_path = config.model_saving_path_topofallfeature if hasattr(config, 'model_saving_path_topofallfeature') else os.path.join(config.root_path_topofallfeature, species, "results/esm/kfold_model_saving/")
    
    # 强制使用ESM模型路径，而不是其他模型路径
    if 'gcn' in model_saving_path.lower() or 'gat' in model_saving_path.lower() or 'sage' in model_saving_path.lower() or 'cnn' in model_saving_path.lower() or 'lstm' in model_saving_path.lower():
        print(f"警告: 配置文件指向非ESM模型路径({model_saving_path})，将替换为ESM模型路径")
        # 从根路径构建到ESM模型的路径
        root_path = config.root_path_topofallfeature
        model_saving_path = os.path.join(root_path, species, "results/esm/kfold_model_saving/")
        print(f"修正后的模型路径: {model_saving_path}")
    
    # 确保模型保存目录存在
    os.makedirs(model_saving_path, exist_ok=True)
    
    # 设置输出目录 - 修改为data/物种/features
    if args.output_dir:
        output_dir = args.output_dir
    else:
        # 获取项目根目录
        project_root = os.path.abspath(os.path.join(current_dir, "../.."))
        # 创建默认输出路径: data/物种/features
        output_dir = os.path.join(project_root, "data", species, "features")
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置参数
    raw_data_path = config.raw_data_path_topofallfeature
    train_batch_size = config.train_batch_size_topofallfeature
    test_batch_size = config.test_batch_size_topofallfeature
    
    print(f"处理物种: {species}, fold: {fold}")
    print(f"模型路径: {model_saving_path}")
    print(f"输出目录: {output_dir}")
    
    # 构建fold路径和模型文件路径
    fold_path = os.path.join(kfold_root_path, f'fold{fold}')
    model_file_name = f'model_dict_for_fold_{fold}.pkl'
    model_file_path = os.path.join(model_saving_path, model_file_name)
    
    if not os.path.exists(model_file_path):
        print(f"错误: 模型文件不存在: {model_file_path}")
        sys.exit(1)
    
    # 加载数据集
    train_data_path = os.path.join(fold_path, 'train_data.txt')
    test_data_path = os.path.join(fold_path, 'test_data.txt')
    
    train_list = pd.read_csv(train_data_path, sep='\t', index_col=False)['Ensembl'].tolist()
    test_list = pd.read_csv(test_data_path, sep='\t', index_col=False)['Ensembl'].tolist()
    
    print(f"训练集大小: {len(train_list)}")
    print(f"测试集大小: {len(test_list)}")
    
    train_data = SeqDataset(gene_list=train_list, raw_data_path=raw_data_path) 
    test_data = SeqDataset(gene_list=test_list, raw_data_path=raw_data_path)  

    # 分割训练集为训练和验证
    train_size = int(0.8 * len(train_data))
    valid_size = len(train_data) - train_size
    train_data, valid_data = torch.utils.data.random_split(train_data, [train_size, valid_size])
    
    # 创建数据加载器
    train_loader = DataLoader(dataset=train_data, batch_size=train_batch_size, shuffle=True, collate_fn=collate_fn)
    valid_loader = DataLoader(dataset=valid_data, batch_size=test_batch_size, shuffle=False, collate_fn=collate_fn)
    test_loader = DataLoader(dataset=test_data, batch_size=test_batch_size, shuffle=False, collate_fn=collate_fn)
    
    # 加载模型
    print(f"加载模型: {model_file_path}")
    try:
        checkpoint = torch.load(model_file_path, map_location=device)
        model = checkpoint
        print("模型结构:")
        print(f"全连接层权重: {model.fc_g2.weight.shape}")
    except Exception as e:
        print(f"加载模型时出错: {e}")
        print("\n推荐解决方案:")
        print("1. 确保模型路径指向ESM模型，而不是其他类型模型")
        print("2. 使用正确的物种名称")
        print("3. 检查配置文件中的model_saving_path_topofallfeature设置")
        sys.exit(1)
    
    # 提取特征
    print("从训练集提取特征...")
    features_array, labels_array = extract_features(model, train_loader, device)
    
    # 保存提取的特征和标签
    feature_file = os.path.join(output_dir, f"{species}_esm_fold{fold}_features.npy")
    label_file = os.path.join(output_dir, f"{species}_esm_fold{fold}_labels.npy")
    np.save(feature_file, features_array)
    np.save(label_file, labels_array)
    print(f"特征已保存到: {feature_file}")
    print(f"标签已保存到: {label_file}")
    
    # 生成可视化图表
    print("生成可视化图表...")
    output_prefix = os.path.join(output_dir, f"{species}_esm_fold{fold}")
    visualize_and_save(features_array, labels_array, output_prefix, species)
    
    print("完成!")

if __name__ == "__main__":
    main()



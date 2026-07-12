#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
特征提取与绘图脚本 - 优化版
支持对多个物种和多个模型进行批量特征提取和绘图
"""

import os
import sys
import time
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# 模型配置 - 可以方便添加或修改支持的模型
MODEL_CONFIGS = {
    "cnn": {
        "script_path": "runners/cnn/cnn_feature.py",
        "max_fold": 4
    },
    "lstm": {
        "script_path": "runners/lstm/lstm_feature.py",
        "max_fold": 4
    },
    "transformer": {
        "script_path": "runners/transformer/transformer_feature.py",
        "max_fold": 4
    },
    "esm": {
        "script_path": "runners/esm/esm_feature.py",
        "max_fold": 2
    },
    "gat": {
        "script_path": "runners/gat/gat_feature.py",
        "max_fold": 2
    }
}

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='批量特征提取与绘图')
    
    # 基本参数
    parser.add_argument('--species', type=str, required=True, nargs='+', 
                        help='要处理的物种列表，例如: celegans fgraminearum')
    parser.add_argument('--models', type=str, nargs='+', 
                        default=["cnn", "lstm", "transformer", "esm", "gat"],
                        help='要运行的模型列表，默认全部')
    
    # 控制参数
    parser.add_argument('--start_fold', type=int, default=0, 
                        help='起始fold索引')
    parser.add_argument('--max_fold', type=int, default=None, 
                        help='最大fold索引，默认使用每个模型的预设值')
    parser.add_argument('--sequential', action='store_true', 
                        help='顺序执行，不使用并行')
    parser.add_argument('--dry_run', action='store_true', 
                        help='只打印命令但不执行')
    parser.add_argument('--python_path', type=str, default='python',
                        help='Python解释器路径')
    
    # 控制台输出
    parser.add_argument('--verbose', action='store_true',
                        help='显示详细输出')
    
    return parser.parse_args()

def generate_commands(args):
    """生成所有需要执行的命令"""
    commands = []
    
    for species in args.species:
        for model_name in args.models:
            if model_name not in MODEL_CONFIGS:
                print(f"警告: 未知模型 {model_name}，跳过")
                continue
                
            model_config = MODEL_CONFIGS[model_name]
            script_path = model_config["script_path"]
            
            # 确定最大fold索引
            max_fold = args.max_fold if args.max_fold is not None else model_config["max_fold"]
            
            # 生成每个fold的命令
            for fold in range(args.start_fold, max_fold + 1):
                cmd = f"{args.python_path} {script_path} --species {species} -fold {fold}"
                
                # 添加元数据用于显示
                commands.append({
                    "command": cmd,
                    "species": species,
                    "model": model_name,
                    "fold": fold
                })
    
    return commands

def execute_command(cmd_info, args):
    """执行单个命令并处理结果"""
    start_time = time.time()
    cmd = cmd_info["command"]
    species = cmd_info["species"]
    model = cmd_info["model"]
    fold = cmd_info["fold"]
    
    # 构建描述
    description = f"{species} - {model} (fold {fold})"
    
    if args.dry_run:
        print(f"[DRY RUN] {description}: {cmd}")
        return True, description, 0
    
    try:
        # 设置环境变量，可能有助于性能优化
        env = os.environ.copy()
        
        # 针对不同模型优化线程设置
        if model == "esm":
            env["OMP_NUM_THREADS"] = "4" 
            env["MKL_NUM_THREADS"] = "4"
        elif model == "transformer":
            env["OMP_NUM_THREADS"] = "2"
            env["MKL_NUM_THREADS"] = "2"
        
        # 执行命令并获取输出
        print(f"开始执行: {description}")
        process = subprocess.run(
            cmd, 
            shell=True, 
            env=env,
            capture_output=not args.verbose,
            text=True,
            check=True
        )
        
        duration = time.time() - start_time
        print(f"✓ 完成: {description} (耗时: {duration:.2f}秒)")
        return True, description, duration
        
    except subprocess.CalledProcessError as e:
        duration = time.time() - start_time
        print(f"✗ 失败: {description} (耗时: {duration:.2f}秒)")
        if args.verbose:
            print(f"错误详情: {e}")
            if e.stdout:
                print(f"标准输出: {e.stdout}")
            if e.stderr:
                print(f"错误输出: {e.stderr}")
        return False, description, duration

def main():
    """主函数"""
    args = parse_arguments()
    
    # 生成所有需要执行的命令
    commands = generate_commands(args)
    
    # 显示任务摘要
    print(f"将处理 {len(args.species)} 个物种, {len(args.models)} 个模型")
    print(f"总共 {len(commands)} 个任务")
    
    # 按物种分组统计
    for species in args.species:
        species_count = sum(1 for cmd in commands if cmd["species"] == species)
        print(f"  - {species}: {species_count} 个任务")
    
    # 按模型分组统计
    for model_name in args.models:
        model_count = sum(1 for cmd in commands if cmd["model"] == model_name)
        print(f"  - {model_name}: {model_count} 个任务")
    
    if args.dry_run:
        print("\n*** 仅打印命令，不执行 ***\n")
    
    # 创建结果记录
    results = {
        "success": [],
        "failed": []
    }
    
    # 开始执行
    start_time = time.time()
    
    # 决定是否并行执行
    if args.sequential:
        print("\n按顺序执行命令...\n")
        for cmd_info in tqdm(commands, desc="总进度"):
            success, description, duration = execute_command(cmd_info, args)
            if success:
                results["success"].append((description, duration))
            else:
                results["failed"].append((description, duration))
    else:
        print("\n并行执行命令...\n")
        # 根据任务类型选择合适的并发数
        # ESM任务内部已经是并行的，所以总体并发较低
        max_workers = 2 if "esm" in args.models else 4
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(execute_command, cmd_info, args) for cmd_info in commands]
            
            for future in tqdm(futures, desc="总进度"):
                success, description, duration = future.result()
                if success:
                    results["success"].append((description, duration))
                else:
                    results["failed"].append((description, duration))
    
    # 计算总执行时间
    total_time = time.time() - start_time
    
    # 打印结果摘要
    print("\n执行结果摘要:")
    print(f"总执行时间: {total_time:.2f} 秒")
    print(f"成功任务: {len(results['success'])}/{len(commands)}")
    print(f"失败任务: {len(results['failed'])}/{len(commands)}")
    
    # 如果有失败任务，打印它们
    if results["failed"]:
        print("\n失败的任务:")
        for desc, duration in results["failed"]:
            print(f"  - {desc} (耗时: {duration:.2f}秒)")
    
    # 按物种和模型分类统计执行时间
    print("\n各物种执行时间:")
    for species in args.species:
        species_time = sum(duration for desc, duration in results["success"] if species in desc)
        print(f"  - {species}: {species_time:.2f}秒")
    
    print("\n各模型执行时间:")
    for model_name in args.models:
        model_time = sum(duration for desc, duration in results["success"] if model_name in desc)
        print(f"  - {model_name}: {model_time:.2f}秒")

if __name__ == "__main__":
    main()

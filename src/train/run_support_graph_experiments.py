import subprocess


def require_epgat_env():
    try:
        import torch  # noqa: F401
        import dgl  # noqa: F401
    except ImportError:
        raise SystemExit("错误：当前环境缺少 torch 或 dgl，请先执行 conda activate EPGAT 再运行。")


def main():
    require_epgat_env()
    print("【support graph 基线实验 runner】")
    print("开始执行训练脚手架...")
    subprocess.check_call(["python", "-m", "src.train.train_support_graph_baseline"])
    print("开始执行评估脚手架...")
    subprocess.check_call(["python", "-m", "src.eval.evaluate_support_graph_baseline"])
    print("support graph baseline 实验 runner 执行完成。")


if __name__ == "__main__":
    main()

# PhaSePred

> 输入蛋白质氨基酸序列，预测其相分离倾向。
> [Chen et al. 2022 *PNAS*](https://doi.org/10.1073/pnas.2115369119) 的开源实现。
> 北京大学李婷婷实验室（Tingting Li Lab）维护。

[English](README.md) · [![ci](https://github.com/NotWhiteBlank/PhaSePred/actions/workflows/ci.yml/badge.svg)](https://github.com/NotWhiteBlank/PhaSePred/actions/workflows/ci.yml)

## PhaSePred 是什么

**相分离（phase separation）** 是蛋白质在细胞内自发凝聚成液滴状区室的过程——像油滴在水中形成。这些"无膜细胞器"（应激颗粒、核仁、Cajal 小体等）无需膜结构就能组织细胞内的生化反应。相分离失调与神经退行性疾病（ALS、阿尔茨海默病）和癌症密切相关。

本工具从氨基酸序列出发，预测蛋白质参与相分离的倾向，并区分两种机制、四个模式：

| 模式 | 含义 | 适用物种 |
|------|------|---------|
| **SaPS** | 自驱动相分离——蛋白自己就能凝聚 | 任意物种（8 个特征） |
| **PdPS** | 伴侣依赖相分离——需要结合伙伴蛋白 | 任意物种（8 个特征） |
| **hSaPS** | 人类自驱动相分离 | 人类蛋白（10 个特征） |
| **hPdPS** | 人类伴侣依赖相分离 | 人类蛋白（10 个特征） |

输出是 0 到 1 之间的分数，分数越高表示越可能属于该模式。

十个特征如下，8 特征模式不含最后两项：

| 特征 | 含义 | 计算来源 | 8f | 10f |
|------|------|---------|:--:|:--:|
| **Hydropathy** | 逐残基归一化 Kyte-Doolittle 疏水均值 | LocalCIDER | ✓ | ✓ |
| **FCR** | 带电残基比例 | LocalCIDER | ✓ | ✓ |
| **IDR** | 内在无序区域占比 | ESpritz（DisProt 模型 `D`、`sw 0`，5% FPR） | ✓ | ✓ |
| **LCR** | 低复杂度区域占比 | SEG | ✓ | ✓ |
| **PScore** | 氨基酸组成复杂度评分 | PScore | ✓ | ✓ |
| **PLAAC** | 类朊病毒结构域 NLLR 分数 | PLAAC | ✓ | ✓ |
| **catGRANULE** | 基于 Bolognesi et al. 2016 的颗粒化打分 | `catgranule` 包 | ✓ | ✓ |
| **DeepCoil** | 卷曲螺旋，按论文取 0.82 阈值二值化 | DeepCoil | ✓ | ✓ |
| **Phos freq** | 磷酸化位点频率 | PhosphoSitePlus | | ✓ |
| **DeepPhase** | DeepPhase 相分离打分 | DeepPhase（随包分发） | | ✓ |

## 快速开始

### 安装

Python 包要求 `>=3.12`，四条路线任选其一：

```bash
pixi install                                                        # pixi
uv venv --python 3.12 && uv pip install .                           # uv
conda env create -f environment.yml && conda activate phasepred     # conda
python -m venv .venv && .venv/bin/pip install -e . -e packages/catgranule   # 原生 venv
```

外部特征工具各自安装：

```bash
bash tools/SEG/install.sh
bash tools/PLAAC/install.sh
bash tools/PScore/install.sh
bash tools/ESpritz/install.sh
bash tools/DeepCoil/install.sh
bash tools/LocalCIDER/install.sh
bash tools/PhosphoSitePlus/install.sh   # 仅 hSaPS / hPdPS 需要
```

然后校验：

```bash
phasepred check-tools          # 表格，每个 MISSING 行给出安装提示
phasepred check-tools --strict # 有必需工具缺失时退出 1，脚本和 CI 用这个
```

`SEG`、`PLAAC`、`LocalCIDER`、`catGRANULE`、`DeepPhase` 随包分发，装完即可用。`PScore`、`ESpritz`、`DeepCoil`、`PhosphoSitePlus` 需要跑上面的安装器。

### 运行

```bash
# 8 特征模式（任意物种）
phasepred predict --fasta proteins.fasta --mode SaPS --output scores.csv
phasepred predict --fasta proteins.fasta --mode PdPS --output scores.csv

# 10 特征模式（人类蛋白）
phasepred predict --fasta proteins.fasta --mode hSaPS --output scores.csv
phasepred predict --fasta proteins.fasta --mode hPdPS --output scores.csv

# 用 UniProt ID 输入
phasepred predict --ids "P35637,Q9Y2W1" --mode SaPS --output scores.csv    # 内联
phasepred predict --ids my_ids.txt --mode SaPS --output scores.csv         # 文件，每行一个
cat my_ids.txt | phasepred predict --ids - --mode SaPS --output scores.csv # stdin

# FASTA 与 ID 可混用
phasepred predict --fasta my.fasta --ids "P35637" --mode SaPS --output scores.csv

# 只算特征，不做预测
phasepred features-from-fasta --input my.fasta --output features.csv
```

输出 CSV 除 `score` 外还包含全部特征列的数值。

## 复现

四张训练表在 `data/processed/`，是论文 Dataset S2/S3 特征列的序列化（base/human × train/test）。

从零训练与 AUC 验证各一条命令，脚本在仓库根目录：

```bash
./train.sh      # 重训 40 个模型，与随包工件逐字节比对
./validate.sh   # 打印四模式 AUC 对论文的对照表
```

`train.sh` 把 40 个模型写进 `runs/retrain/`，再与 `src/phasepred/data/models/` 逐个做 sha256 比对，最后打印 `N/40 byte-identical`。退出码 `0` 全部相同、`1` 有工件不同、`2` 预检失败。逐字节复现的条件是 `xgboost>=3.2,<3.3`。

`validate.sh` 跑三层验证（论文层 AUC、泄露、产物一致性），把对照表打印到屏幕，报告写进 `runs/validation/`：

```
MODE        PAPER      MEASURED      DELTA  VERDICT
SaPS-8      0.862      0.860729    -0.0013  pass
PdPS-8      0.739      0.737884    -0.0011  pass
hSaPS-10    0.924      0.917515    -0.0065  pass
hPdPS-10    0.827      0.816171    -0.0108  fail
```

退出码 `0` 无硬门禁失败、`1` 硬门禁失败、`2` 输入或用法错误。指标与划分落在 `products/A_paper_split_recomputed/`。

## 引用

> Chen, Z., Hou, C., Wang, L., Yu, C., Chen, T., Shen, B., Hou, Y.,
> Li, P., Li, T. (2022). Screening membraneless organelle participants
> with machine-learning models that integrate multimodal features.
> *Proceedings of the National Academy of Sciences* 119(24), e2115369119.
> https://doi.org/10.1073/pnas.2115369119

catGRANULE 打分实现源自 Bolognesi et al. 2016（*Cell Reports* 16:222-231）。

## License

本仓库（代码、模型、文档）采用 **MIT License**，详见 [LICENSE](LICENSE)。

> Copyright (c) 2026 Tingting Li Lab, Department of Biochemistry and Molecular Biology,
> School of Basic Medical Sciences, Peking University.

第三方组件及其许可见 [NOTICE](NOTICE)。

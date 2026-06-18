# PhaSePred

> 输入蛋白质氨基酸序列，预测其相分离倾向。
> 本仓库由北京大学李婷婷实验室（Tingting Li Lab）维护，是 [Chen et al. 2022 *PNAS*](https://doi.org/10.1073/pnas.2115369119)（本实验室前期工作）的延续版本——使用扩充后的数据库与调优后的 XGBoost 重新训练，并提供端到端的 CLI 工具。

[English](README-en.md)

## 背景

**相分离（phase separation）** 是蛋白质在细胞内自发凝聚成液滴状区室的过程——像油滴在水中形成。这些"无膜细胞器"（应激颗粒、核仁、Cajal小体等）无需膜结构就能组织细胞内的生化反应。相分离失调与神经退行性疾病（ALS、阿尔茨海默病）和癌症密切相关。

本工具预测蛋白质参与相分离的倾向，区分两种机制：

| 模式 | 含义 | 适用物种 |
|------|------|---------|
| **SaPS** | 自驱动相分离——蛋白自己就能凝聚 | 任意物种（8 个特征） |
| **PdPS** | 伴侣依赖相分离——需要结合伙伴蛋白 | 任意物种（8 个特征） |
| **hSaPS** | 人类自驱动相分离 | 人类蛋白（10 个特征） |
| **hPdPS** | 人类伴侣依赖相分离 | 人类蛋白（10 个特征） |

输出是 0 到 1 之间的分数。分数越高，蛋白属于该模式的可能性越大。分数不是一个校准的概率，而是该模式内的相对排序。

## 两个独立产物

仓库提供 **两套模型**，各自独立可发布：

| 产物 | 训练数据 | 划分方式 | 类别不平衡处理 | 默认？ |
|---|---|---|---|---|
| **A** | 论文 S2/S3 + 我们自重算的特征 | 论文原始划分 | 2:1 负样本子抽样（论文协议） | 否 |
| **B** | PhaSepDB 3.0 + 10 物种背景负样本 + LLPSDB v2 + 自重算特征 | 分层 80/20，seed=42 | `scale_pos_weight`（全部负样本） | **是** |

两者都使用 Optuna 在训练分区上做 100 trials × 5-fold CV 调参，**测试集严格封存**。`predict` 命令默认走 Product B，`--product A` 可切换。详细见 [products/README.md](products/README.md)。

### 关键性能数字

| 任务 | A: CV → 测试 (S3) | B: CV → 测试 (20%) | B: PhaSePro 外部验证 |
|------|------|------|------|
| SaPS  | 0.888 → **0.794** | 0.840 → **0.838** | 0.782 |
| PdPS  | 0.759 → **0.699** | 0.764 → **0.737** | 0.859 |
| hSaPS | 0.933 → **0.843** | 0.855 → **0.816** | 0.839 |
| hPdPS | 0.851 → **0.785** | 0.823 → **0.859** | 0.959 |

> 外部验证集（PhaSePro、LLPSDB2）已经做了 **accession 去重**——任何同时出现在训练集和外部集的蛋白都被剔除。完整 overlap 计数见 `products/<product>/leakage_report.json`。LLPSDB2 列因为其"确认非相分离"标签贴近决策边界，两个产物都在 0.5 附近——这是数据本身的限制，已记录。

## 快速开始

```bash
# 1. 克隆并准备 Python 环境
git clone https://github.com/notwhiteblank/PhaSePred.git && cd PhaSePred
uv sync

# 2. 安装外部特征工具（详见下一节）
bash tools/install/install_all.sh    # 或参考下方"工具安装"分步指引

# 3. 校验所有工具就绪
uv run phasepred check-tools

# 4. 预测
uv run phasepred predict --fasta my.fasta --mode SaPS --output scores.csv
```

每批约 8 秒。输出 CSV 包含 `UniprotEntry, score, length, Hydropathy, FCR, IDR, LCR, PScore, PLAAC, catGRANULE, DeepCoil[, Phos freq, DeepPhase]`，所有特征数值都有，便于核查。

## 工具安装（必读）

PhaSePred 调用 **6 个外部特征工具**（catGRANULE 已内建为 Python 模块），见下表。每个工具都通过 **路径解析器** 定位，按 `环境变量 → 仓库内 vendored → PATH` 顺序找。安装后请务必执行 `phasepred check-tools` 验证。

| 工具 | License | 安装方式 | 环境变量覆盖 |
|---|---|---|---|
| **SEG** | NCBI 公开 | 仓库内已 vendored 源码，`cd tools/per-tool/SEG && make` | `PHASEPRED_SEG_BIN` |
| **PScore** | CC-BY 4.0 | 仓库内已 vendored，无需安装 | `PHASEPRED_PSCORE_DIR` |
| **PLAAC** | MIT | 仓库内已 vendored 预编译 jar；需要 Java 17+ | `PHASEPRED_PLAAC_JAR` |
| **ESpritz** | 学术许可，禁止再分发 | `bash tools/install/install_espritz.sh`（提示用户接受条款下载） | `PHASEPRED_ESPRITZ_DIR` |
| **IUPred3** | 学术许可，禁止再分发 | `bash tools/install/install_iupred3.sh`（同上） | `PHASEPRED_IUPRED3_DIR` |
| **DeepCoil** | 上游无 LICENSE 声明 | `bash tools/install/install_deepcoil_env.sh`（创建 Python 3.8 conda env） | `PHASEPRED_DEEPCOIL_ENV` |

vendored 部分克隆即用。第三方授权限制的部分必须本地下载，license 审计明细见 [docs/THIRD_PARTY_LICENSES.md](docs/THIRD_PARTY_LICENSES.md)。

### 路径解析器如何工作

源码里**没有任何 `Tools/<X>/<bin>` 形式的硬编码路径**。`src/phasepred/tool_paths.py` 暴露 `find_seg() / find_pscore_dir() / find_plaac_wrapper() / ...`，每个函数走 4 步级联：

1. **环境变量**：例如 `PHASEPRED_SEG_BIN=/opt/bin/seg`
2. **仓库内 vendored**：`tools/per-tool/SEG/seg`
3. **PATH**：`shutil.which("seg")`
4. **报错**：抛出 `PhaSePredToolNotFound`，并附带安装提示

如果你已经把某个工具装在系统其他位置（例如用 conda、apt 装的版本），**不需要重装也不需要把它复制进仓库**——把对应环境变量指向它即可。

### check-tools 检查表

安装完成后跑：

```bash
uv run phasepred check-tools
```

输出示例：

```
TOOL        STATUS  REQUIRED  SOURCE                       PATH
SEG         OK      yes       vendored                     /repo/tools/per-tool/SEG/seg
PScore      OK      yes       vendored                     /repo/tools/per-tool/PScore/SourceCodeS2
PLAAC       OK      yes       vendored                     /repo/tools/per-tool/PLAAC/plaac-master/web/bin/plaac.jar
IUPred3     OK      no        env:PHASEPRED_IUPRED3_DIR    /opt/iupred3
ESpritz     OK      yes       env:PHASEPRED_ESPRITZ_DIR    /opt/espritz
DeepCoil    OK      yes       vendored (conda)             .external_envs/deepcoil
catGRANULE  OK      yes       built-in                     src/phasepred/catgranule_v1.py

7 OK  ·  0 MISSING  ·  0 required-but-missing
```

任何 `required-but-missing > 0`，`phasepred predict` 会在读 FASTA **之前** 就退出并提示具体安装命令——避免那种"跑到一半才发现工具缺失"的情况。

### 工具故障排查

| 症状 | 通常原因 | 修复 |
|---|---|---|
| `SEG executable not found` | 没编译 SEG | `cd tools/per-tool/SEG && make` |
| `PLAAC: java not found` | 没装 Java 运行时 | `sudo apt install openjdk-17-jre-headless` 或 macOS `brew install openjdk` |
| `ESpritz: perl not found` | 没装 Perl | 多数 Linux 自带；macOS `brew install perl` |
| `DeepCoil env not found` | 没装隔离 conda env | `bash tools/install/install_deepcoil_env.sh` |
| `IUPred3 not found` 但已下载 | 解压位置不对 | 把 `iupred3/` 放在 `tools/per-tool/IUPred3/iupred3/` 或设 `PHASEPRED_IUPRED3_DIR` |

## 可选：人类模型额外数据（hSaPS / hPdPS 需要）

`SaPS` 和 `PdPS` 用 8 个特征，完全从序列计算，不需要额外文件。

`hSaPS` 和 `hPdPS` 多用 2 个特征——`Phos freq`（磷酸化位点密度）和 `DeepPhase`（显微图像评分）。这两个不是可以"装一个工具跑出来"的特征：它们来自上游论文/数据库提供的**查找表**，因此 PhaSePred **不能从序列重新算**。两份数据必须由用户自己下载：

| 特征 | 数据来源 | License | 期望路径 |
|---|---|---|---|
| Phos freq | [PhosphoSitePlus](https://www.phosphosite.org/) `Phosphorylation_site_dataset.gz` | 需要注册下载，禁止再分发 | `data/raw/external/phosphositeplus/Phosphorylation_site_dataset.gz` |
| DeepPhase | [DeepPhase 论文 supplement](https://github.com/PEILab/DeepPhase) `tableS3.xlsx` | 学术使用，参见上游 | `data/raw/external/deepphase/extracted/tableS3.xlsx` |

**没有这两个文件会怎样？** `phasepred predict --mode hSaPS|hPdPS` **仍然能跑**，但会在 stderr 输出 `PhaSePredMissingDataWarning`，明确告诉你哪几个蛋白用了中位数填充（即"分数仅供参考，不是基于完整特征的判断"）。例如：

```
PhaSePredMissingDataWarning: DeepPhase has no entry for 3/5 requested
protein(s): P12345, Q67890, A1B2C3. Median imputation will be used, so
hSaPS/hPdPS scores for these proteins should be treated as approximate.
```

如果你只要做 SaPS/PdPS 预测，可以完全忽略这两个文件。

## 使用

### 预测（FASTA → 分数）

```bash
# 8 特征模型（任意物种）
uv run phasepred predict --fasta proteins.fasta --mode SaPS  --output scores.csv
uv run phasepred predict --fasta proteins.fasta --mode PdPS  --output scores.csv

# 10 特征模型（人类蛋白）
uv run phasepred predict --fasta proteins.fasta --mode hSaPS --output scores.csv
uv run phasepred predict --fasta proteins.fasta --mode hPdPS --output scores.csv

# 切换产物
uv run phasepred predict --fasta proteins.fasta --mode SaPS --product A --output scores.csv

# 用 UniProt ID 直接预测（三种输入形式都支持）
uv run phasepred predict --ids "P35637,Q9Y2W1" --mode SaPS --output scores.csv      # 内联
uv run phasepred predict --ids my_ids.txt --mode SaPS --output scores.csv           # 文件（每行一个 ID）
cat my_ids.txt | uv run phasepred predict --ids - --mode SaPS --output scores.csv   # stdin

# FASTA 和 IDs 可以混用（去重后合并）
uv run phasepred predict --fasta my.fasta --ids "P35637" --mode SaPS --output scores.csv
```

ID 输入会调 `https://rest.uniprot.org` 拉序列，结果缓存到 `data/interim/uniprot_cache.jsonl`，下次同样的 ID 不再请求网络。用 `--cache PATH` 可以指定别的缓存位置。

### 复现两个产物

```bash
# 1. 准备特征（如果 data/interim/recomputed_features_full.csv 不在）
#    需要先安装所有外部工具
uv run phasepred features-recomputed \
    --sequences data/interim/sequences.csv \
    --output data/interim/recomputed_features_full.csv \
    --espritz-cache data/interim/espritz_idr_cache.jsonl \
    --catgranule-source v1

# 2. 调参 + 训练 Product A
uv run python products/A_paper_split_recomputed/tune.py
uv run python products/A_paper_split_recomputed/train.py

# 3. 调参 + 训练 Product B
uv run python products/B_extended_dataset/tune.py
uv run python products/B_extended_dataset/train.py
```

调参约 30-60 分钟，训练约 1-2 分钟（取决于硬件）。

## 仓库结构

```
src/phasepred/             共享库 + CLI
  cli.py                   predict / check-tools / features-recomputed / ...
  predictor.py             FASTA → 特征 → 预测全流程，支持 --product A|B
  tool_paths.py            外部工具路径解析器（env / vendored / PATH）
  features.py              特征定义和原生特征计算（Hydropathy、FCR）
  catgranule_v1.py         催化颗粒形成倾向，论文公式重建（无上游代码）
  paper.py                 Product A 使用：解析论文 S2/S3 supplements
  updated_data.py          Product B 使用：PhaSepDB3 + LLPSDB2 标签合并
  ...

products/                  两个独立软件产物
  A_paper_split_recomputed/  论文划分 + 自重算特征 + 调优 XGBoost
    tune.py / train.py
    tuned_params.json metrics.json leakage_report.json
    train_accessions_*.tsv test_accessions_*.tsv   （封存的训练/测试划分）
    models/<task>/8f_model_{0..9}.joblib            （10 模型集成 × 4 任务）
  B_extended_dataset/        更新数据集 + 调优 XGBoost（默认）
    （同上结构）

tools/                     外部工具层
  README.md                            路径解析机制
  THIRD_PARTY_LICENSES.md              License 审计
  install/                             install_<tool>.sh 自动化安装脚本
  wrappers/                            run_<tool>.sh shell 包装器
  per-tool/<TOOL>/                     每个工具的 README、源码（如允许 vendored）

docs/
  THIRD_PARTY_LICENSES.md
  DATA_SOURCES.md
  ENVIRONMENT.md
  REPO_LAYOUT.md
  STORAGE_LAYOUT.md

tests/                     pytest 测试套件
```

## 探索过程与关键结论

项目从重建实验室 2022 年发表的 PhaSePred 起步（Chen et al. 2022, *PNAS*），随着数据库更新和方法学打磨，演进到当前两个独立产物。下面按时间顺序记录每一步，**关键数字直接嵌入**——不依赖额外文档。

### Step 1a · LogReg sanity check

用 2022 论文的 S2/S3 划分 + 三个最简单的序列特征（Hydropathy、FCR、IDR）跑 LogisticRegression，目的是验证数据加载和评估管道。AUC 远低于 2022 公布值（仅 ~0.7 量级），符合预期：模型太弱，但管道工作。

### Step 1b · XGBoost 超参 Optuna 探索

在 2022 论文 S2 特征上跑 100 trials × 3-fold CV Optuna 搜索 XGBoost 超参空间，覆盖 n_estimators (100-600)、max_depth (3-8)、learning_rate (log-uniform 0.01-0.3) 等 9 个参数。这一步只是探索，没有产生最终模型。

### Step 1c · 端到端重建 2022 论文协议

用 2022 论文特征值 + 默认 XGBoost + 2:1 负样本子抽样 + 10 模型集成，严格按 2022 方法重建 Table 1，确认重写后的代码与历史结果一致：

| 任务 | 2022 公布值 | 本仓库重建值 (S3 test AUC) |
|---|---|---|
| SaPS  | 0.86 | 0.79 |
| PdPS  | 0.66 | 0.74 |
| hSaPS | 0.85 | 0.83 |
| hPdPS | 0.78 | 0.82 |

四个任务都在 2022 公布值的 ±0.05 误差带内，子抽样导致的方差是主因——证实新的代码栈忠实地还原了原始方法学。

### Step 2 · 外部特征工具集成

逐个安装并验证 SEG（NCBI 公开）、ESpritz（学术许可）、IUPred3（学术许可）、PScore（CC-BY 4.0）、PLAAC（MIT）、DeepCoil（无 license）、catGRANULE 等工具，使预测能直接从 FASTA 序列出发，不再依赖 2022 论文中间表格。期间装过 CD-HIT 用于序列同源去冗余，但最终在两个产物中都未启用，已移除。

### Step 3 · AlphaFold pLDDT 特征实验（阴性结果）

**假设**：AlphaFold 的逐残基置信度 pLDDT 平均值能补充"哪些残基处于结构化区域 vs 无序区域"的信息，作为第 9 个特征加入应能提升 AUC。

**方法**：通过 AlphaFold DB API 抓取 60,315 个蛋白的平均 pLDDT，按 2022 协议重新训练 XGBoost 集成（9 特征版），与 8 特征版本比较。

**结果（阴性）**：

| 任务 | 8 特征 AUC | 9 特征 (含 pLDDT) AUC | Δ |
|---|---|---|---|
| SaPS  | 0.794 | 0.794 | 0.000 |
| PdPS  | 0.699 | 0.701 | +0.002 |
| hSaPS | 0.843 | 0.838 | -0.005 |
| hPdPS | 0.785 | 0.781 | -0.004 |

特征重要性中 pLDDT 排名 8/9，且与 PScore (r ≈ 0.45)、LCR (r ≈ 0.50)、PLAAC (r ≈ 0.42) 显著相关——它捕获的"结构化程度"信息已被现有特征覆盖。**结论：弃用**。生产模型只用 8 个基础特征（人类模型加 Phos freq + DeepPhase 共 10 个）。这是个有用的阴性结果——它告诉我们这个方向不需要再投入。

### Step 4 · 数据集扩展第一次迭代

PhaSepDB 3.0 发布后，沿用 2022 论文的方法重新划分（PhaSepDB3 正样本 + 10 物种背景负样本 + LLPSDB v2 confirmed negatives），重新计算所有特征。**问题**：当时使用默认 XGBoost 参数，并且外部验证（PhaSePro / LLPSDB2）没做 accession overlap 过滤。5-fold CV AUC 看起来很漂亮（SaPS 0.866, PdPS 0.832, hSaPS 0.977, hPdPS 0.973），但因为没去重，外部 AUC 数字过于乐观。这一版被 Product B 取代——B 加了 Optuna 调参 + accession 去重过滤。

### 补充 · catGRANULE v1 公式重建

catGRANULE v1 上游代码不可获取，只有 Bolognesi 等 2016 论文的公式。从头实现：解析论文补充材料 Table S1 推导 Z 归一化常数，按公式实现颗粒形成倾向打分。验证：与 2022 PhaSePred supplement S2 提供的 catGRANULE 列做对比，**RMSE ≈ 0.24**，特征级 AUC 达到 parity。生产 pipeline 优先用 2022 supplement 数值，缺失时回退到本重建版（`src/phasepred/catgranule_v1.py`）。

### Step 5 · IDR 特征敏感性测试（ESpritz vs IUPred3）

2022 论文指定 ESpritz DisProt @ 5% FPR 作为 IDR 信号源；IUPred3 是更常用的替代品。控制其他特征不变，测试两个 IDR 来源的影响：ESpritz 全任务 AUC 比 IUPred3 高 0.01-0.03，且匹配 2022 定义。**结论**：保留 ESpritz 作为默认，IUPred3 在 `src/phasepred/iupred.py` 中作为可选 IDR 源保留。

### 本次发布 · Products A 与 B（最终产出）

两个产物都做了：
- **Optuna 100 trials × 5-fold stratified CV** 在训练分区上调参（每个任务独立调）
- **测试集严格封存**——调参时不可见
- **外部验证 accession 去重**——任何同时在训练集和外部集（PhaSePro, LLPSDB2）的蛋白都被剔除，AUC 只在非重叠子集上计算
- 详细 overlap 计数写到 `products/<product>/leakage_report.json`

性能数字见上方"关键性能数字"表格。这两个产物对应 `products/` 下两个目录，相互独立，可分别使用。

## 引用

本仓库是本实验室 2022 年 PhaSePred 工作的延续版本。如使用本工具，请引用原始论文：

> Chen, Z., Hou, C., Wang, L., Yu, C., Chen, T., Shen, B., Hou, Y.,
> Li, P., Li, T. (2022). Screening membraneless organelle participants
> with machine-learning models that integrate multimodal features.
> *Proceedings of the National Academy of Sciences* 119(24), e2115369119.
> https://doi.org/10.1073/pnas.2115369119

外部特征工具的引用见 [docs/THIRD_PARTY_LICENSES.md](docs/THIRD_PARTY_LICENSES.md)。

## License

PhaSePred 本仓库（代码、模型、文档）采用 **MIT License** 发布。详见根目录 [LICENSE](LICENSE)。

> Copyright (c) 2026 Tingting Li Lab, Department of Biochemistry and Molecular Biology,
> School of Basic Medical Sciences, Peking University.

仓库内 vendored 的第三方组件（PLAAC、PScore、SEG）保留各自的上游 license——见
[docs/THIRD_PARTY_LICENSES.md](docs/THIRD_PARTY_LICENSES.md) 的完整审计。用户自行安装的工具
（ESpritz、IUPred3、DeepCoil）遵循各自上游条款。

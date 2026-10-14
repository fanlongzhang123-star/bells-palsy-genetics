# 代码与稿件内容的对应关系 (MANIFEST)

论文 Code Availability 承诺的内容与本包脚本的对应：

| 论文中的分析 | 脚本 | 产出 |
|---|---|---|
| 协调与 meta 分析 | `scripts/meta_analysis.py` | `out/meta_bp_r13_plus_gcst90018843.csv`；正文 OR 1.217 / beta +0.1960 |
| SuSiE 可信集处理 | `scripts/fine_mapping_susie.py` | 三个可信集，CS1 四变异 / CS2 RUNX2 |
| MHC 内显著变异分布 | `scripts/mhc_snp_scan.py` | 表 1；图 4a |
| 187 个经典 HLA 等位基因 | `scripts/hla_alleles.py` | 表 2；图 4b |
| 共定位 | `scripts/colocalisation_and_credible_set.py` | 表 3；PP.H4 及自校验 PP.H4 = 1.000 |
| 可信集位点直接查询 | `scripts/colocalisation_and_credible_set.py` | 图 3；四个变异在 11 个表型上的效应 |
| 跨性状孟德尔随机化 | `scripts/cross_trait_mr.py` | 补充图 S1b；正文 OR 与 95% CI |
| 全表型组扫描 | `scripts/phenome_scan.py` | 表 4；图 4c |
| **UK Biobank（Pan-UKB ICD-10 G51）复制查证** | `scripts/ukb_replication.py` | 正文表 2；补充表 S2；`out/ukb_G51_credibleset_lookup.csv`（输入：`results/credibleset_variants_finngen_r13.tsv`） |
| z 分数相关 | `scripts/genetic_relatedness.py` | 补充图 S1a |
| 图 1 数据 | `scripts/extract_figure_data.py` | `work/fig_manhattan.tsv`, `work/fig_region.tsv` |
| 全部图件数据汇总 | `scripts/build_figure_data.py` | `figures/data/*.tsv` |
| 全部图件绘制 | `figures.R` | Fig. 1-4、Fig. S1、Graphical abstract |

## 说明

- `results/` 内是上述脚本产出的小体积结果表，便于第三方在不重跑全流程的情况下核对数值。
- `coloclib.py`、`susie.py`、`mrlib.py` 是为本环境无法安装 scipy/pandas 而自行实现的统计分布与估计量，未使用第三方统计包。
- 计算环境：Python 3.14（仅标准库 + numpy）、R 4.6.1（base graphics，无外部包）。
- 稿件排版与校验脚本（qa_manuscript.py 等）不属于科学分析，未包含在本包内，可向通讯作者索取。
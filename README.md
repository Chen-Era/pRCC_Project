# pRCC_Project

面向肾癌（pRCC）的单细胞与空间转录组一体化分析脚本集合，覆盖从原始矩阵导入、质控、标准化/降维、聚类与标记物识别、批次整合、功能富集、通路分析、CNV 推断、轨迹推断、TF 调控网络（SCENIC）到基因集富集（GSVA/QuSAGE）等关键环节。

本 README 基于仓库内现有脚本的真实参数与输出更新，提供依赖、CLI 用法与最小工作示例，帮助快速复现与组合使用各模块。

## 目录与脚本概览

- Seruat3Normalization.R：单细胞数据导入、过滤、CellCycle 评分、归一化、可变基因、Scale、PCA/UMAP/t-SNE，保存处理后的 Seurat 对象（RDS）。
- Seurat3Cluster.R：基于 Seurat 的聚类、簇统计与标记基因识别，输出多个汇总表与 RDS。
- MNN.R：使用 batchelor::fastMNN 进行批次整合，并进行 UMAP/t-SNE、聚类与保存结果。
- ClusterProfile-GO.R：使用 clusterProfiler 进行 GO 富集（BP/MF/CC），输出 CSV。
- ClusterProfile-Pathway.R：使用 clusterProfiler 进行 KEGG 通路富集，输出 CSV。
- GSVA.ssGSEA.R：对表达矩阵进行 ssGSEA（GSVA 包），输出打分矩阵。
- Monocle2.Pseudotime.R：基于 Monocle2 的伪时间轨迹分析，从 Seurat RDS 生成轨迹图。
- inferCNV.R：基于已计算的 CNV 矩阵与分组信息，判定肿瘤细胞并输出细胞/簇层面的 CNV 指标与图形。
- QuSAGE.R：对 Seurat RDS 的各聚类进行 QuSAGE 基因集活性分析，输出密度曲线与置信区间图及表格。
- SCENIC.sh：基于 pySCENIC 的 GRN 推断、motif 富集与 AUCell 步骤封装脚本。

## 环境与依赖（按脚本）

通用建议：R 4.x；建议使用 Bioconductor 安装生信包；Python 环境安装 pySCENIC 及其依赖。部分脚本仅检测依赖而不自动安装，如提示“Missing packages”请先手动安装。

- Seruat3Normalization.R：Seurat、dplyr、ggplot2 等（随 Seurat 生态）。
- Seurat3Cluster.R：SingleCellExperiment、scater、plyr、reshape2、Seurat、mclust、dplyr。
- MNN.R：batchelor、SingleCellExperiment、scater、Seurat（以及绘图所需依赖）。
- ClusterProfile-GO.R：clusterProfiler、org.Hs.eg.db、dplyr（需要 Bioconductor）。
- ClusterProfile-Pathway.R：clusterProfiler、org.Hs.eg.db、dplyr（需要 Bioconductor）。
- GSVA.ssGSEA.R：GSVA、limma、stringr。
- Monocle2.Pseudotime.R：monocle、Seurat、reshape2、dplyr、ggplot2（保存图）。
- inferCNV.R：ggplot2。
- QuSAGE.R：qusage、Seurat。
- SCENIC.sh：pySCENIC（arboreto、ctx 等），需准备对应物种的 motif/cistarget 数据库文件。

R 依赖安装示例（部分包需 Bioconductor）：

- install.packages(c("Seurat", "dplyr", "ggplot2", "plyr", "reshape2", "mclust", "limma", "stringr", "ggplot2"))
- if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
- BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GSVA", "qusage", "monocle", "SingleCellExperiment", "scater", "batchelor"))

Python/pySCENIC 安装示例：

- pip install pyscenic

## CLI 用法与最小示例

以下示例默认输出目录不存在时会自动创建，且大多脚本将结果写入指定 output_dir 或 ./Result。

1) 归一化与降维（Seruat3Normalization.R）
- 用法：Rscript Seruat3Normalization.R <MatrixDir> <sample.txt> <output_dir> <prefix> <genecolumn> <cellMinGene> <cellMaxGene> <geneExpMinCell> <MaxMTPercent> <regressOut>
- 说明：
  - MatrixDir：10x Matrix 或稀疏矩阵目录
  - sample.txt：样本元数据
  - prefix：输出前缀
  - genecolumn：基因列名
  - 其余为过滤与回归参数
- 输出：处理后的 Seurat 对象 RDS（含 PCA/UMAP/t-SNE）。

2) 聚类与标记物（Seurat3Cluster.R）
- 用法：Rscript Seurat3Cluster.R <input_rds> <prefix> <output_dir> <resolution> <marker_method> <logFC_threshold> <min.pct> <only.pos>
- 输出：
  - <prefix>_GraphClust.Statistics.txt
  - <prefix>_GraphClust.Summary_Cell.txt
  - <prefix>_GraphClust.AllMarkerGenes.txt
  - <prefix>_GraphClust.seuset.rds

3) 批次整合（MNN.R）
- 用法：Rscript MNN.R <input_rds> <output_dir> <minCell> <k> <resolution>
- 输出：整合后的对象与降维结果（UMAP/t-SNE）及聚类结果（具体文件名以脚本输出为准）。

4) GO 富集（ClusterProfile-GO.R）
- 用法：Rscript ClusterProfile-GO.R <gene_list.txt> [output_dir]
- 输出：GO_enrichment.csv（包含 BP/MF/CC 三类）

5) KEGG 富集（ClusterProfile-Pathway.R）
- 用法：Rscript ClusterProfile-Pathway.R <gene_list.txt> [output_dir] [organism:hsa]
- 输出：KEGG_enrichment.csv

6) ssGSEA（GSVA.ssGSEA.R）
- 用法：Rscript GSVA.ssGSEA.R <gene_set.gmt|txt> <expression_matrix.txt> [output_dir]
- 输出：ssgsea.txt（行为样本，列为基因集得分）

7) 伪时间（Monocle2.Pseudotime.R）
- 用法：Rscript Monocle2.Pseudotime.R <seuset.rds> <output_dir> [gene_file]
- 输出：trajectory_state.png、trajectory_cluster.png、trajectory_pseudotime.png

8) CNV 指标（inferCNV.R）
- 用法：Rscript inferCNV.R <observations.txt> <references.txt> <group.txt> <output_dir>
- 输出：
  - CNV_cells.txt、CNV_clusters.txt
  - correlation_against_meanSquare.pdf、CNVPlot.png

9) QuSAGE 基因集活性（QuSAGE.R）
- 用法：Rscript QuSAGE.R <input_rds> <geneset_files_comma_separated> <output_dir>
- 说明：geneset 文件可为 GMT 或两列 txt（type, gene），多个以逗号分隔
- 输出：在 <output_dir>/QuSAGE/<基因集名>/ 下生成每个簇的 table_*.txt、plotDC_*.png/pdf、plotCI_*.png

10) SCENIC（SCENIC.sh, 基于 pySCENIC）
- 典型流程包含三步：基因调控网络推断（grn）、motif 富集（ctx）与 AUCell 打分（aucell）。
- 建议运行：bash SCENIC.sh -h 查看完整参数与数据库路径配置；需准备对应物种的 feather 排序数据库与 motif 注释。
- 输出：SCENIC 相关中间与最终结果文件将保存在指定输出目录（包括候选调控网络、富集结果与细胞打分矩阵等）。

## 输入/输出格式约定

- 基因列表：一列 SYMBOL 文本；GO/KEGG 会自动映射到 ENTREZID（物种为人：org.Hs.eg.db）。
- 表达矩阵：制表符分隔文本（基因为行或列视脚本而定，GSVA/ssGSEA 要求列名为样本/细胞）。
- Seurat 对象：RDS 格式；推荐使用本仓库提供的归一化脚本生成，保证槽位兼容。
- CNV 矩阵：空格分隔文本；observations 为肿瘤候选，references 为参考细胞；group.txt 两列（细胞ID, 分组/簇）。

## 注意事项与最佳实践

- 依赖安装：部分脚本仅检测、不会自动安装依赖；如报 Missing packages，请先按上文安装。
- 版本兼容：Seurat v3/v4 槽位命名存在差异；QuSAGE/Monocle2 已做基本兼容处理，但仍建议使用本仓库流程生成的对象。
- 性能与资源：SCENIC、MNN 等计算量较大，建议先在小样本子集验证；必要时增加内存与并行线程。
- 路径与编码：路径含空格请加引号；确保输入文件为 UTF-8，无 BOM。
- 可复现性：建议在脚本外部设定 set.seed() 并记录会话信息；输出目录建议按阶段维护 ./Result/ 与 ./Test/（先小规模验证再全量运行）。

## 参考运行顺序（示例）

1) 归一化/降维 → 2) 聚类/标记物 → 3) 批次整合（如需） → 4) 功能富集/通路 → 5) GSVA/QuSAGE → 6) 轨迹/SCENIC → 7) CNV 指标整合与可视化。

根据课题需要可灵活调整与迭代，每个模块均可独立运行与组合。

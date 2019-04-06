## ----eval = TRUE, message=FALSE------------------------------------------
library(MetaNeighbor)
library(SummarizedExperiment)
data(mn_data)
data(GOmouse)

## ----eval=TRUE,fig.width=4,fig.height=3, results='hide'------------------
AUROC_scores = MetaNeighbor(dat = mn_data,
                            experiment_labels = as.numeric(factor(mn_data$study_id)),
                            celltype_labels = metadata(colData(mn_data))[["cell_labels"]],
                            genesets = GOmouse,
                            bplot = TRUE)

## ----eval= TRUE----------------------------------------------------------
head(AUROC_scores)

## ----eval = TRUE---------------------------------------------------------
library(MetaNeighbor)
data(mn_data)

## ----eval = TRUE---------------------------------------------------------
var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
head(var_genes)

## ----eval = TRUE---------------------------------------------------------
length(var_genes)

## ----eval=TRUE-----------------------------------------------------------
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = mn_data,
                             study_id = mn_data$study_id,
                             cell_type = mn_data$cell_type)

## ----eval=TRUE,fig.width=7,fig.height=6.5--------------------------------
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)

## ----eval = TRUE---------------------------------------------------------
top_hits = topHits(cell_NV = celltype_NV,
                   dat = mn_data,
                   study_id = mn_data$study_id,
                   cell_type = mn_data$cell_type,
                   threshold = 0.9)
top_hits

## ----eval = TRUE, message=FALSE------------------------------------------
library(MetaNeighbor)
library(SummarizedExperiment)
data(mn_data)
data(GOmouse)

## ----eval = TRUE,fig.width=4,fig.height=3, results='hide'----------------
AUROC_scores = MetaNeighbor(dat = mn_data,
                            experiment_labels = as.numeric(factor(mn_data$study_id)),
                            celltype_labels = metadata(colData(mn_data))[["cell_labels"]],
                            genesets = GOmouse,
                            bplot = TRUE,
                            fast_version = TRUE)

## ----eval = TRUE, fig.width = 7, fig.height = 6.5------------------------
var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = mn_data,
                             study_id = mn_data$study_id,
                             cell_type = mn_data$cell_type,
                             fast_version = TRUE)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)

## ----eval = TRUE, message = FALSE----------------------------------------
library(SingleCellExperiment)
library(Matrix)
baron <- readRDS('baron-human.rds')
segerstolpe <- readRDS('segerstolpe.rds')

## ----eval = TRUE---------------------------------------------------------
common_genes <- intersect(rownames(baron), rownames(segerstolpe))
baron <- baron[common_genes,]
segerstolpe <- segerstolpe[common_genes, !(segerstolpe$cell_type1 %in% c('not applicable', 'co-expression'))]

## ----eval = TRUE---------------------------------------------------------
new_colData = data.frame(
  study_id = rep(c('baron', 'segerstolpe'), c(ncol(baron), ncol(segerstolpe))),
  cell_type = c(as.character(colData(baron)$cell_type1), colData(segerstolpe)$cell_type1)
)
pancreas <- SingleCellExperiment(
  Matrix(cbind(assay(baron, 1), assay(segerstolpe, 1)), sparse = TRUE),
  colData = new_colData
)
dim(pancreas)
rm(baron); rm(segerstolpe)

## ----eval = TRUE, fig.width=7,fig.height=6.5-----------------------------
var_genes = variableGenes(dat = pancreas, exp_labels = pancreas$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = pancreas,
                             study_id = pancreas$study_id,
                             cell_type = pancreas$cell_type,
                             fast_version = TRUE)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)

## ----eval = TRUE---------------------------------------------------------
all_pancreas <- readRDS('all_pancreas.rds')
dim(all_pancreas)

## ----eval = TRUE, eval=TRUE,fig.width=7,fig.height=6.5-------------------
var_genes = variableGenes(dat = all_pancreas, exp_labels = all_pancreas$Study_ID)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = all_pancreas,
                             study_id = all_pancreas$Study_ID,
                             cell_type = all_pancreas$Celltype,
                             fast_version = TRUE)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)

## ----eval = TRUE,fig.width=4,fig.height=3, results='hide'----------------
data(GOhuman)
small_pancreas = all_pancreas[, all_pancreas$Celltype %in% c('alpha', 'beta', 'delta')]
celltype_matrix = model.matrix(~small_pancreas$Celltype - 1)
colnames(celltype_matrix) = levels(as.factor(small_pancreas$Celltype))
AUROC_scores = MetaNeighbor(dat = small_pancreas,
                            experiment_labels = as.numeric(factor(small_pancreas$Study_ID)),
                            celltype_labels = celltype_matrix,
                            genesets = GOhuman,
                            bplot = TRUE,
                            fast_version = TRUE)


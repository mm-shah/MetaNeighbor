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


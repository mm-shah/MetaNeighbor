#' @name MetaNeighbor_sample_data
#' @title MetaNeighbor_sample_data
#'
#' @description A list containing a collection of following: a gene matrix,
#' cell type labels, experiment labels and sets of genes.
#'
#' @docType data
#'
#' @format
#' \describe{
#'   \item{data}{A gene-by-sample expression matrix consisting of 3157 rows
#'   (genes) and 1051 columns (cell types)}
#'   \item{exp.lab}{A numerical vector of length 1051 that indicates the source
#'   of each sample (1 = Zeisel et al, 2 = Tasic et al)}
#'   \item{cell.lab}{1051x1 binary matrix that indicates whether a cell belongs
#'   to the SstNos cell type (1=yes, 0 = no)}
#'   \item{genesets}{List containing gene symbols for 10 GO function (GO:0016853
#'   , GO:0005615, GO:0005768, GO:0007067, GO:0065003, GO:0042592, GO:0005929,
#'   GO:0008565, GO:0016829, GO:0022857) downloaded from the Gene Ontology
#'   Consortuum August 2015 \url{http://www.geneontology.org/}}
#' }
#'
#' @source Dataset:\url{https://github.com/mm-shah/MetaNeighbor/tree/master/data} 1. Zeisal et al. \url{http://science.sciencemag.org/content/347/6226/1138}
#' 2. Tasic et al. \url{http://www.nature.com/neuro/journal/v19/n2/full/nn.4216.html}
#' @rdname MetaNeighbor_sample_data

"MetaNeighbor_sample_data"

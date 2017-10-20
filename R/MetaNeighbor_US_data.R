#' @name MetaNeighbor_US_data
#' @title MetaNeighbor_US_data
#'
#' @description A dataset containing a collection of following: a gene matrix,
#' cell type labels, experiment labels and sets of genes.
#'
#' @docType data
#'
#' @format
#' \describe{
#'   \item{data}{A gene-by-sample expression matrix consisting of 3157 rows
#'   (genes) and 1051 columns (samples)}
#'   \item{pheno}{A sample metadata table (1051x2), that lists the dataset and
#'   cell type for each sample with column names: Study_ID and Celltype}
#'   \item{celltypes}{A vector of length 39 that lists the names of all cell
#'   types}
#'   }
#'
#' @source Dataset:\url{https://github.com/mm-shah/MetaNeighbor/tree/master/data}
#'
#' @rdname MetaNeighbor_US_data

"MetaNeighbor_US_data"

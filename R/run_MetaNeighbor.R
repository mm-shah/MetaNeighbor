#' Runs MetaNeighbor
#'
#' For each gene set of interest, the function builds a network of rank
#' correlations between all cells. Next,It builds a network of rank correlations
#' between all cells for a gene set. Next, the neighbor voting predictor
#' produces a weighted matrix of predicted labels by performing matrix
#' multiplication between the network and the binary vector indicating cell type
#' membership, then dividing each element by the null predictor (i.e., node
#' degree). That is, each cell is given a score equal to the fraction of its
#' neighbors (including itself), which are part of a given cell type. For
#' cross-validation, we permute through all possible combinations of
#' leave-one-dataset-out cross-validation, and we report how well we can recover
#' cells of the same type as area under the receiver operator characteristic
#' curve (AUROC). This is repeated for all folds of cross-validation, and the
#' mean AUROC across folds is reported. Calls
#' \code{\link{neighbor_voting_LeaveOneExpOut}}.
#'
#' @param data A gene-by-sample expression matrix.
#' @param experiment_labels A numerical vector that indicates the source of each
#' sample.
#' @param celltype_labels A matrix that indicates the cell type of each sample.
#' @param genesets Gene sets of interest provided as a list of vectors.
#' @param file_ext Specify the desired path and file name for output file.
#' @return There are three outputs of the method:
#'      i.) A vector of AUROC scores representing the mean for each gene set
#'      tested is returned directly
#'      (see \code{\link{neighbor_voting_LeaveOneExpOut}}).
#'     ii.) A list containing the AUROC score for each cell type across all
#'     folds of cross-validation can be found in the first output file.
#'     iii.) A matrix containing the means for each cell type across folds of
#'     cross-validaiton can be found in the second output file.
#'
#' @examples
#' data(MetaNeighbor_sample_data)
#' AUROC_scores = run_MetaNeighbor(data = MetaNeighbor_sample_data$data,
#'                                 experiment_labels = MetaNeighbor_sample_data$exp.lab,
#'                                 celltype_labels = MetaNeighbor_sample_data$cell.lab,
#'                                 genesets = MetaNeighbor_sample_data$genesets,
#'                                 file_ext = "filename")
#' hist(AUROC_scores,
#'      main = "Sst Chodl",
#'      xlab = "AUROC Scores",
#'      breaks = 10,
#'      xlim = c(0,1))
#' abline(v = mean(AUROC_scores), col = "red", lty = 2, lwd = 2)
#'
#' @export
#'

run_MetaNeighbor <-function(data,
                            experiment_labels,
                            celltype_labels,
                            genesets,
                            file_ext) {

    ROCs              <- vector(mode = "list", length = length(genesets))
    names(ROCs)       <- names(genesets)
    nv_mat            <- matrix(data = 0,
                                ncol = dim(celltype_labels)[2],
                                nrow = length(genesets))
    rownames(nv_mat)  <- names(genesets)
    colnames(nv_mat)  <- colnames(celltype_labels)

    for (l in 1:length(genesets)){
        print(l)
        geneset     <- genesets[[l]]
        m           <- match(rownames(data), geneset)
        dat_sub     <- data[!is.na(m),]
        dat_sub     <- stats::cor(dat_sub, method = "s")
        dat_sub     <- as.matrix(dat_sub)
        rank_dat    <- dat_sub
        rank_dat[]  <- rank(dat_sub, ties.method = "average", na.last = "keep")
        rank_dat[is.na(rank_dat)] <- 0
        rank_dat    <- rank_dat/max(rank_dat)
        ROCs[[l]]   <- neighbor_voting_LeaveOneExpOut(experiment_labels,
                                                        celltype_labels,
                                                        rank_dat,
                                                        means = FALSE)
    }

    for(i in 1:length(ROCs)){
        nv_mat[i,] <- rowMeans(ROCs[[i]][[1]], na.rm = TRUE)
    }

    save(ROCs, file=paste(file_ext, "IDscore.list.Rdata", sep="."))
    save(nv_mat, file=paste(file_ext, "IDscore.matrix.Rdata", sep="."))
    return(rowMeans(nv_mat, na.rm = TRUE))
}

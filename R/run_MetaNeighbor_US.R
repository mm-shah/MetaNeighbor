#' Runs unsupervised version of MetaNeighbor
#'
#' When it is difficult to know how cell type labels compare across datasets this
#' function helps users to make an educated guess about the overlaps without
#' requiring in-depth knowledge of marker genes
#'
#' @param vargenes A vector of high variance genes.
#' @param data A gene-by-sample expression matrix.
#' @param pheno A sample metadata table, that lists the dataset and cell type
#' for each sample with column names: Study_ID and Celltype.
#'
#' @return The output is a cell type-by-cell type mean AUROC matrix, which is built by
#' treating each pair of cell types as testing and training data for MetaNeighbor, then
#' taking the average AUROC for each pair (NB scores will not be identical because each
#' test cell type is scored out of its own dataset, and the differential heterogeneity
#' of datasets will influence scores).
#'
#' @examples
#' data(MetaNeighbor_US_data)
#' var_genes = get_variable_genes(MetaNeighbor_US_data$data, MetaNeighbor_US_data$pheno)
#' celltype_NV = run_MetaNeighbor_US(var_genes,
#'                                   MetaNeighbor_US_data$data,
#'                                   MetaNeighbor_US_data$pheno)
#' celltype_NV
#'
#' @export
#'

run_MetaNeighbor_US <- function(vargenes, data, pheno){

    pheno$StudyID_CT <- paste(pheno$Study_ID, pheno$Celltype, sep = "+")
    celltypes <- unique(pheno$StudyID_CT)
    cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
    rownames(cell_labels) <-colnames(data)
    colnames(cell_labels) <- celltypes

    for(i in 1:length(celltypes)){
        type <- celltypes[i]
        matching_celltype <- match(pheno$StudyID_CT, type)
        cell_labels[!is.na(matching_celltype),i]  <- 1
    }

    matching_vargenes <- match(rownames(data), vargenes)
    matching_vargenes_count   <- length(which(is.na(matching_vargenes)== FALSE))

    if(matching_vargenes_count < 2){
        stop("matching_vargenes should have more than 1 matching genes!", call. = TRUE)
    } else if(matching_vargenes_count < 5) {
        warning("matching_vargenes should have more matching genes!", immediate. = TRUE)
    }

    cor_data    <- stats::cor(data[!is.na(matching_vargenes),], method="s")
    rank_data   <- cor_data*0
    rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
    rank_data[is.na(rank_data)] <- 0
    rank_data   <- rank_data/max(rank_data)
    sum_in      <- (rank_data) %*% cell_labels
    sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum), ncol = dim(sum_in)[2], nrow = dim(sum_in)[1])
    predicts    <- sum_in/sum_all

    cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
    colnames(cell_NV) <- colnames(cell_labels)
    rownames(cell_NV) <- colnames(cell_labels)

    for(i in 1:dim(cell_labels)[2]){
        predicts_temp <- predicts

        matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
        unique_studyID    <- unique(pheno[!is.na(matching_celltype), "Study_ID"])
        matching_studyID  <- match(pheno$Study_ID, unique_studyID)
        pheno2            <- pheno[!is.na(matching_studyID),]
        predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
        predicts_temp     <- apply(abs(predicts_temp), MARGIN = 2, FUN = rank, na.last= "keep", ties.method="average")


        filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
        matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
        filter[!is.na(matches),1:length(celltypes)] <- 1

        negatives = which(filter == 0, arr.ind = TRUE)
        positives = which(filter == 1, arr.ind = TRUE)

        predicts_temp[negatives] <- 0

        np <- colSums(filter, na.rm = TRUE)
        nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=TRUE))
        p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm = TRUE)

        cell_NV[i,]= (p/np - (np+1)/2)/nn
    }

    cell_NV <- (cell_NV+t(cell_NV))/2
    return(cell_NV)
}

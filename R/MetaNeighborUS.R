#' Runs unsupervised version of MetaNeighbor
#'
#' When it is difficult to know how cell type labels compare across datasets this
#' function helps users to make an educated guess about the overlaps without
#' requiring in-depth knowledge of marker genes
#'
#' @param var_genes vector of high variance genes.
#' @param data SummarizedExperiment object containing gene-by-sample
#' expression matrix. 
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#' 
#' @return The output is a cell type-by-cell type mean AUROC matrix, which is 
#' built by treating each pair of cell types as testing and training data for 
#' MetaNeighbor, then taking the average AUROC for each pair (NB scores will not
#' be identical because each test cell type is scored out of its own dataset, 
#' and the differential heterogeneity of datasets will influence scores).
#'
#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' celltype_NV = MetaNeighborUS(var_genes = var_genes, 
#'                              dat = mn_data, 
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' celltype_NV
#'
#' @export
#'

MetaNeighborUS <- function(var_genes, dat, i = 1, study_id, cell_type){
    
    dat    <- SummarizedExperiment::assay(dat, i = i)
    samples <- colnames(dat)
    
    #check obj contains study_id
    if(length(study_id)!=length(samples)){
        stop('study_id length does not match number of samples')
    }
    
    #check obj contains cell_type
    if(length(cell_type)!=length(samples)){
        stop('cell_type length does not match number of samples')
    } 
    
    pheno <- as.data.frame(cbind(study_id,cell_type), stringsAsFactors = FALSE)
    pheno$StudyID_CT <- paste(pheno$study_id, pheno$cell_type, sep = "|")
    celltypes   <- unique(pheno$StudyID_CT)
    cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
    rownames(cell_labels) <-colnames(data)
    colnames(cell_labels) <- celltypes

    for(i in seq_along(celltypes)){
        type <- celltypes[i]
        matching_celltype <- match(pheno$StudyID_CT, type)
        cell_labels[!is.na(matching_celltype),i]  <- 1
    }

    matching_vargenes <- match(rownames(dat), var_genes)
    matching_vargenes_count   <- sum(!is.na(matching_vargenes))

    if(matching_vargenes_count < 2){
        stop("matching_vargenes should have more than 1 matching genes!",
             call. = TRUE)
    } else if(matching_vargenes_count < 5) {
        warning("matching_vargenes should have more matching genes!", 
                immediate. = TRUE)
    }

    cor_data    <- stats::cor(dat[!is.na(matching_vargenes),], method="s")
    rank_data   <- cor_data*0
    rank_data[] <- rank(cor_data, ties.method = "average", na.last = "keep")
    rank_data[is.na(rank_data)] <- 0
    rank_data   <- rank_data/max(rank_data)
    sum_in      <- (rank_data) %*% cell_labels
    sum_all     <- matrix(apply(rank_data, MARGIN = 2, FUN = sum), 
                          ncol = dim(sum_in)[2], 
                          nrow = dim(sum_in)[1])
    predicts    <- sum_in/sum_all

    cell_NV     <- matrix(0, ncol=length(celltypes), nrow=length(celltypes))
    colnames(cell_NV) <- colnames(cell_labels)
    rownames(cell_NV) <- colnames(cell_labels)

    for(i in seq_len(dim(cell_labels)[2])){
        predicts_temp <- predicts

        matching_celltype <- match(pheno$StudyID_CT, colnames(cell_labels)[i])
        unique_studyID    <- unique(pheno[!is.na(matching_celltype),"study_id"])
        matching_studyID  <- match(pheno$study_id, unique_studyID)
        pheno2            <- pheno[!is.na(matching_studyID),]
        predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
        predicts_temp     <- apply(abs(predicts_temp), 
                                   MARGIN = 2, 
                                   FUN = rank, 
                                   na.last= "keep", 
                                   ties.method="average")


        filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
        matches <- match(pheno2$StudyID_CT, colnames(cell_labels)[i])
        filter[!is.na(matches),seq_along(celltypes)] <- 1

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

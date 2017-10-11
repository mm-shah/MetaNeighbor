run_MetaNeighbor_US <- function(vargenes, data, celltypes, pheno){

  cell_labels <- matrix(0, ncol=length(celltypes), nrow=dim(pheno)[1])
  rownames(cell_labels) <-colnames(data)
  colnames(cell_labels) <- celltypes

  for(i in 1:length(celltypes)){
    type <- celltypes[i]
    matching_celltype <- match(pheno$Celltype, type)
    cell_labels[!is.na(matching_celltype),i]  <- 1
  }

  matching_vargenes <- match(rownames(data), vargenes)
  cor_data    <- cor(data[!is.na(matching_vargenes),], method="s")
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

    matching_celltype <- match(pheno$Celltype, colnames(cell_labels)[i])
    unique_studyID    <- unique(pheno[!is.na(matching_celltype), "Study_ID"])
    matching_studyID  <- match(pheno$Study_ID, unique_studyID)
    pheno2            <- pheno[!is.na(matching_studyID),]
    predicts_temp     <- predicts_temp[!is.na(matching_studyID),]
    predicts_temp     <- apply(abs(predicts_temp), MARGIN = 2, FUN = rank, na.last= "keep", ties.method="average")

    filter  <- matrix(0, ncol=length(celltypes), nrow=dim(pheno2)[1])
    matches <- match(pheno2$Celltype, colnames(cell_labels)[i])
    filter[!is.na(matches),1:length(celltypes)] <- 1

    negatives = which(filter == 0, arr.ind = T)
    positives = which(filter == 1, arr.ind = T)

    predicts_temp[negatives] <- 0

    np <- colSums(filter, na.rm=T)
    nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x==0,na.rm=T))
    p  <- apply(predicts_temp, MARGIN = 2, FUN = sum, na.rm=T)

    cell_NV[i,]= (p/np - (np+1)/2)/nn
  }

  cell_NV <- (cell_NV+t(cell_NV))/2
  return(cell_NV)

}

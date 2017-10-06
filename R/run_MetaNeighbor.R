run_MetaNeighbor <- function(data, experiment_labels, celltype_labels, genesets, file_ext) {

  ROCs              <- vector(mode = "list", length = length(genesets))
  names(ROCs)       <- names(genesets)
  nv_mat            <- matrix(data = 0, ncol = dim(celltype_labels)[2], nrow = length(genesets))
  rownames(nv_mat)  <- names(genesets)
  colnames(nv_mat)  <- colnames(celltype_labels)

  for (l in 1:length(genesets)){
    print(l)
    geneset     <- genesets[[l]]
    m           <- match(rownames(data), geneset)
    dat_sub     <- data[!is.na(m),]
    dat_sub     <- cor(dat_sub, method = "s")
    dat_sub     <- as.matrix(dat_sub)
    rank_dat    <- dat_sub
    rank_dat[]  <- rank(dat_sub, ties.method = "average", na.last = "keep")
    rank_dat[is.na(rank_dat)] <- 0
    rank_dat    <- rank_dat/max(rank_dat)
    ROCs[[l]]   <- neighbor_voting_LeaveOneExpOut(experiment_labels, celltype_labels, rank_dat, means = F)
  }

  for(i in 1:length(ROCs)){
    nv_mat[i,] <- rowMeans(ROCs[[i]][[1]], na.rm = T)
  }

  save(ROCs, file=paste(file_ext, "IDscore.list.Rdata", sep="."))
  save(nv_mat, file=paste(file_ext, "IDscore.matrix.Rdata", sep="."))
  return(rowMeans(nv_mat, na.rm = T))
}


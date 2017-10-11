get_top_hits <- function(cell_NV, pheno, threshold=0.95, filename) {

  type_by_study <- table(pheno[ , c("Celltype", "Study_ID")])
  m <- match(rownames(cell_NV), rownames(type_by_study))
  f_a <- !is.na(m)
  f_b <- m[f_a]
  cell_NV <- cell_NV[f_a,f_a]
  type_by_study <- type_by_study[f_b,]

  for(i in 1:dim(type_by_study)[2]){
    filt <- type_by_study[,i] != 0
    cell_NV[filt,filt] <- 0
  }

  diag(cell_NV) <-0
  temp <- vector()

  for(i in 1:dim(cell_NV)[1]){
    temp <- c(temp, which.max(cell_NV[i,]))
  }

  temp <- cbind(rownames(cell_NV), temp)
  for(i in 1:dim(cell_NV)[1]){
    temp[i,2]=cell_NV[i,as.numeric(temp[i,2])]
  }

  recip <- temp[duplicated(temp[,2]),]
  filt  <- as.numeric(temp[,2]) >= threshold
  recip <- rbind(recip,temp[filt,])
  recip <- cbind(recip,c(rep("Reciprocal_top_hit",each=dim(recip)[1]-sum(filt)),rep(paste("Above",threshold,sep="_"),each=sum(filt))))
  recip <- recip[!duplicated(recip[,2]),]

  recip2  <- cbind(rownames(recip),recip[,1:3])
  colnames(recip2) <- c("Celltype_1","Celltype_2","Mean_AUROC","Match_type")
  rownames(recip2) <- NULL

  recip     <- recip2[order(recip2[,3],decreasing=T),]
  recip2    <- as.data.frame(recip)
  recip2[,3]<- round(as.numeric(as.character(recip2[,3])),2)
  write.table(recip, file = filename, sep="\t", quote = F)
  return(recip2)
}

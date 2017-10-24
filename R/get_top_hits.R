#' Find reciprocal top hits
#'
#' Identifies reciprocal top hits and high scoring cell type pairs
#'
#' @param cell_NV A matrix of celltype-to-celltype AUROC scores
#' (output from \code{\link{run_MetaNeighbor_US}})
#' @param pheno A sample metadata table, that lists the dataset and cell type
#' for each sample with column names "Study_ID" and "Celltype".
#' @param threshold Default value 0.95. Must be between [0-1]
#' @param filename Set path and name of output file
#'
#' @return Creates a file and returns a dataframe with cell types that are
#' either reciprocal best matches, and/or those with AUROC values greater than
#' or equal to threshold value

#' @examples
#' data(MetaNeighbor_US_data)
#' var_genes = get_variable_genes(MetaNeighbor_US_data$data,
#'                                MetaNeighbor_US_data$pheno)
#' celltype_NV = run_MetaNeighbor_US(var_genes,
#'                                   MetaNeighbor_US_data$data,
#'                                   MetaNeighbor_US_data$pheno)
#' top_hits = get_top_hits(celltype_NV,
#'                         MetaNeighbor_US_data$pheno,
#'                         threshold=0.9,
#'                         filename="filename.txt")
#' top_hits
#'
#' @export
#'

get_top_hits <- function(cell_NV, pheno, threshold=0.95, filename) {

    pheno$StudyID_CT <- paste(pheno$Study_ID, pheno$Celltype, sep = "+")
    type_by_study <- table(pheno[,"StudyID_CT"])
    m <- match(rownames(cell_NV), rownames(type_by_study))
    f_a <- !is.na(m)
    f_b <- m[f_a]
    cell_NV <- cell_NV[f_a,f_a]
    type_by_study <- type_by_study[f_b]

    # remove within-dataset scores
    for(i in unique(pheno$Study_ID)){
        filt <- grepl(i, row.names(type_by_study != 0))
        cell_NV[filt,filt] <- 0
    }

    # remove self-scores
    diag(cell_NV) <-0
    temp <- vector()

    # identify top hits
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
    recip <- cbind(recip, c(rep("Reciprocal_top_hit",
                                each=dim(recip)[1]-sum(filt)),
                            rep(paste("Above",threshold,sep="_"),
                                each=sum(filt))))
    recip <- recip[!duplicated(recip[,2]),]

    #tidy results
    recip2  <- cbind(rownames(recip),recip[,1:3])
    colnames(recip2) <- c("Celltype_1","Celltype_2","Mean_AUROC","Match_type")
    rownames(recip2) <- NULL

    recip     <- recip2[order(recip2[,3],decreasing=TRUE),]
    recip2    <- as.data.frame(recip)
    recip2[,3]<- round(as.numeric(as.character(recip2[,3])),2)
    utils::write.table(recip, file = filename, sep="\t", quote = FALSE)
    return(recip2)
}

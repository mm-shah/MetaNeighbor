#' Find reciprocal top hits
#'
#' Identifies reciprocal top hits and high scoring cell type pairs
#'
#' @param cell_NV matrix of celltype-to-celltype AUROC scores
#' (output from \code{\link{MetaNeighborUS}})
#' @param dat a SummarizedExperiment object containing gene-by-sample
#' expression matrix.
#' @param i default value 1; non-zero index value of assay containing the matrix
#' data
#' @param study_id a vector that lists the Study (dataset) ID for each sample
#' @param cell_type a vector that lists the cell type of each sample
#' @param threshold default value 0.95. Must be between [0,1]
#'
#' @return Function returns a dataframe with cell types that are either reciprocal best 
#' matches, and/or those with AUROC values greater than or equal to threshold 
#' value

#' @examples
#' data(mn_data)
#' var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
#' celltype_NV = MetaNeighborUS(var_genes = var_genes, 
#'                              dat = mn_data, 
#'                              study_id = mn_data$study_id,
#'                              cell_type = mn_data$cell_type)
#' top_hits = topHits(cell_NV = celltype_NV,
#'                    dat = mn_data,
#'                    study_id = mn_data$study_id,
#'                    cell_type = mn_data$cell_type,
#'                    threshold = 0.9)
#' top_hits
#' 
#' @import gplots RColorBrewer
#' @export
#'

topHits <- function(cell_NV, dat, i = 1, study_id, cell_type, threshold=0.95){
    
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
    
    ct=rownames(cell_NV)
    ct_study = paste(study_id,cell_type,sep='|')
    pheno <- as.data.frame(cbind(study_id,cell_type,ct_study), stringsAsFactors = FALSE)
    type_by_study=table(pheno[,c("ct_study","study_id")])
    m<-match(ct,rownames(type_by_study))
    f.a=!is.na(m)
    f.b=m[f.a]
    cell_NV=cell_NV[f.a,f.a]
    type_by_study=type_by_study[f.b,]
    
    for(i in seq_along(1:dim(type_by_study)[2])){
        filt=type_by_study[,i]!=0
        cell_NV[filt,filt]=0
    }
    
    diag(cell_NV)=0
    
    hits=list()
    names(hits)=
        for(i in seq_along(1:dim(cell_NV)[1])){
            hits[[i]]=which(cell_NV[i,]==max(cell_NV[i,]))
        }
    
    temp=matrix(0,ncol=4,nrow=length(hits))
    for(i in seq_along(1:length(hits))){
        if(!is.na(temp[i,1])){
            a = hits[[i]]
            if(hits[[a]]==i){
                temp[i,]=c(ct[i],ct[a],round(cell_NV[i,a],2),"Reciprocal_top_hit")
                temp[a,] = rep(NA,4)
            }
            else if(cell_NV[i,a]>=threshold){
                temp[i,]=c(ct[i],ct[a],round(cell_NV[i,a],2),paste("Above_threshold:",threshold,sep=''))
            }
            else {
                temp[i,] = rep(NA,4)
            }
        } else {}
    }
    temp = temp[!is.na(temp[,1]),]
    colnames(temp)=c("Study_ID|Celltype_1","Study_ID|Celltype_2","Mean_AUROC"," Match_type")
    ord=order(temp[,3])
    temp=temp[rev(ord),]
    rownames(temp)=c(1:dim(temp)[1])
    temp=as.data.frame(temp)
    return(temp)
}

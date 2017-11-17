#' Evaluates SummarizedExperiment Object provided as input
#' 
#' This script runs checks on the object to verify presense of required data at 
#' its specified location in the object. If any of the checks fail, the program 
#' is terminated with error message. If all checks pass then the control is
#' transfered back to the calling statement.
#' 
#' @param obj A summarizedExperiment object
#' 
#' @return It does not return anything but halts the execution of program if any
#' of the check statements fail.
#' 
#' @examples  
#' data(mn_data)
#' eval_obj(mn_data)
#' 
#' @export
#'

eval_obj <- function(obj){
    samples <- colnames(obj)
    
    #check obj contains sample_id
    if(length(obj$sample_id)!=0){
        if(length(obj$sample_id)!=length(samples)){
            stop('sample_id length does not match number of samples')
        } 
    }else {
        stop(paste(deparse(substitute(obj)),'should contain vector named sample_id'))
    }
    
    #check obj contains study_id
    if(length(obj$study_id)!=0){
        if(length(obj$study_id)!=length(samples)){
            stop('study_id length does not match number of samples')
        } 
    }else {
        stop(paste(deparse(substitute(obj)),'should contain vector named study_id'))
    }
    
    #check obj contains more than 1 unique study_id
    if(length(unique(obj$study_id)) < 2){
        stop('Found only 1 unique study_id. Please use data from more than 1 study!')
    }
    
    #check obj contains cell_type
    if(length(obj$cell_type)!=0){
        if(length(obj$cell_type)!=length(samples)){
            stop('cell_type length does not match number of samples')
        } 
    }else {
        stop(paste(deparse(substitute(obj)),'should contain vector named cell_type'))
    }
    
    #check obj contains cell_labels
    if(length(obj@colData@metadata$cell_labels) == 0){
        stop('cell_labels matrix not present under @metadata section of @colData')
    }

    #check obj contains genesets
    if(length(obj@metadata$genesets) == 0){
        stop('genesets not present under @metadata section')
    }

    #check genesets matches more than 1 genes in gene_matrix
    genes_in_geneset <- as.character(unlist(obj@metadata$genesets))
    genes_in_matrix <- rownames(obj)
    if(length(intersect(genes_in_geneset,genes_in_matrix))<1)
        stop('No matching genes between geneset and gene_matrix')
}
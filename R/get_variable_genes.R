#' Identify a highly variable gene set
#'
#' Identifies genes with high variance compared to their median expression
#' (top quartile) within each experimentCertain function
#'
#' @param data A gene-by-sample expression matrix.
#' @param pheno A sample metadata table, that lists the dataset and cell type
#' for each sample with column names "Study_ID" and "Celltype".
#'
#' @return The output is a vector of gene names that are highly variable in
#' every experiment (intersect)
#'
#' @examples
#' data(MetaNeighbor_US_data)
#' var_genes = get_variable_genes(MetaNeighbor_US_data$data,
#'                                MetaNeighbor_US_data$pheno)
#' var_genes
#'
#' @export
#'

get_variable_genes <- function(data, pheno) {
    var_genes1 <- vector("list")
    experiment <- unique(pheno$Study_ID)
    j <- 1

    for(exp in experiment){
        data_subset   <- data[ , pheno$Study_ID == exp]
        genes_list    <- vector("list")
        median_data   <- apply(data_subset, MARGIN = 1, FUN = stats::median)
        variance_data <- apply(data_subset, MARGIN = 1, FUN = stats::var)
        quant_med     <- unique(stats::quantile(median_data,
                                                probs = seq(from = 0,
                                                            to = 1,
                                                            length = 11),
                                                type = 5))
        genes_list    <- vector("list", length = length(quant_med))

        for(i in 1:length(quant_med)){
            if(i == 1){
            filt1     <- median_data <= quant_med[i]
            var_temp  <- variance_data[filt1]
            quant_var <- stats::quantile(var_temp, na.rm = TRUE)
            filt2     <- var_temp > quant_var[4]
            genes_list[[i]] <- names(var_temp)[filt2]
            } else {
            filt1     <- median_data <= quant_med[i] & median_data > quant_med[i-1]
            var_temp  <- variance_data[filt1]
            quant_var <- stats::quantile(var_temp, na.rm = TRUE)
            filt2     <- var_temp > quant_var[4]
            genes_list[[i]] <- names(var_temp)[filt2]
            }
        }

        temp <- length(genes_list)
        var_genes1[[j]] <- unlist(genes_list[1:temp-1])
        j <- j+1
    }
    var_genes <- Reduce(intersect, var_genes1)
    return(var_genes)
}

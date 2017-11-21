## ----eval = FALSE--------------------------------------------------------
#  library(SummarizedExperiment)
#  exp_mat <- gene_by_sample_matrix
#  genesets <- list_of_gene_sets
#  cell_labels <- sample_by_celltype_matrix
#  colData <- DataFrame(exp_label = vector_of_experiment_labels_of_samples,
#                       sample_id = vector_of_sample_ids_of_samples,
#                       study_id = vector_of_study_id_of_samples,
#                       cell_type = vector_of_cell_type_of_samples,
#                       row.names = colnames(exp_mat)
#                       )
#  data <- SummarizedExperiment(assays=list(gene_matrix=exp_mat),
#                               colData=colData
#                               metadata = list(genesets = genesets))
#  data@colData@metadata <- list(cell_labels = cell_labels)


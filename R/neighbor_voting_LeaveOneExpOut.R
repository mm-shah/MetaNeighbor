neighbor_voting_LeaveOneExpOut <- function( exp_labels, cell_labels, network, means=T){

  # genes.label : needs to be in 1s and 0s
  l <- dim(cell_labels)[2]
  c <- dim(cell_labels)[1]
  e <- unique(exp_labels)


  #print("Make genes label CV matrix")
  test_cell_labels <- matrix(cell_labels, nrow=c, ncol = length(e)*l)
  exp_cols <- rep(e, each = l)

  for (i in 1:length(e)){
    d <- which(exp_labels == i)
    a <- which(exp_cols == i)
    test_cell_labels[d,a] <- 0
  }

  #print("Get sums - mat. mul.")
  #sumin    = ( t(network) %*% test.genes.labels)
  sum_in  <- (network %*% test_cell_labels)

  #print("Get sums - calc sumall")
  sum_all <- matrix(apply(network, MARGIN = 2, FUN = sum), ncol = dim(sum_in)[2], nrow = dim(sum_in)[1])

  #print("Get sums - calc predicts")
  predicts <- sum_in/sum_all

  #print("Hide training data")
  nans <- which(test_cell_labels == 1, arr.ind = T)
  predicts[nans] <- NA

  #Hide other experiment data
  for (i in 1:length(e)){
    d <- which(exp_labels != i)
    a <- which(exp_cols == i)
    predicts[d,a] <- NA
  }

  #print("Rank test data")
  predicts <- apply(abs(predicts), MARGIN = 2, FUN = rank, na.last = "keep", ties.method = "average")
  filter <- matrix(cell_labels, nrow = c, ncol = length(e)*l)

  for (i in 1:length(e)){
    d<-which(exp_labels != i)
    a<-which(exp_cols ==i )
    filter[d,a] <- NA
  }

  negatives <- which(filter == 0, arr.ind = T)
  positives <- which(filter == 1, arr.ind = T)
  predicts[negatives] <- 0

  #print("Calculate ROC - np")
  np <- colSums(filter, na.rm = T) # Postives

  #print("Calculate ROC - nn")
  nn <- apply(filter, MARGIN = 2, FUN = function(x) sum(x == 0, na.rm = T))     # Negatives

  #print("Calculate ROC - p")
  p  <- apply(predicts, MARGIN = 2, FUN = sum, na.rm = T)

  #print("Calculate ROC - rocN")
  rocNV           <- (p/np - (np+1)/2)/nn
  rocNV           <- matrix(rocNV, ncol = length(e), nrow = l)
  colnames(rocNV) <- e
  rownames(rocNV) <- colnames(cell_labels)

  if(means==T){
    scores = list(rowMeans(rocNV ,na.rm = T))
  } else {
    scores = list(rocNV)
  }
}


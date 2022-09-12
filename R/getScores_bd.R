getScores_bd <- function(i, dat, As){
  # scores <- dat[[i]]%*%As[[i]]
  scores <- blockmult(dat[[i]],As[[i]], onmemory = T)
  colnames(scores[[1]]) <- paste0("comp", 1:ncol(scores[[1]]))
  scores[[1]]
}

getScores_bd <- function(i, dat, As){
  # scores <- dat[[i]]%*%As[[i]]
  scores <- blockmult(dat[[i]],As[[i]], onmemory = T)
  colnames(scores) <- paste0("comp", 1:ncol(scores))
  scores
}

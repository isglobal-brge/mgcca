getScores <- function(i, dat, As){
  scores <- dat[[i]]%*%As[[i]]
  colnames(scores) <- paste0("comp", 1:ncol(scores))
  scores
}

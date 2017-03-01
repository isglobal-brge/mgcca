getScores <- function(i, dat, As){
  scores <- dat[[i]]%*%As[[i]]
  scores
}

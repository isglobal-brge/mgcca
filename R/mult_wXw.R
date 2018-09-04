mult_wXw <- function(X, w){
  # X matrix
  # w diag of a matrix
  wX <- mult_wX(X,w)
  wXw <- mult_Xw(wX,w)
  wXw
}

mult_wX <- function(X, w){
  # X matrix
  # w diag of a matrix
  n <- ncol(X)
  v <- rep(w, n)
  wX <- X*v
  wX
}

mult_Xw <- function(X, w){
  # X matrix
  # w diag of a matrix
  n <- nrow(X)
  v <- rep(w, each=n)
  Xw <- X*v
  Xw
}




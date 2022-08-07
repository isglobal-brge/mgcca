mult_wXw_bd <- function(X, w){
  # X matrix
  # w diag of a matrix
  wX <- BigDataStatMeth::bdwproduct( X, w, "wX")
  wXw <- BigDataStatMeth::bdwproduct( wX, w, "Xw")

  return(wXw)
}


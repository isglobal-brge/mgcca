n <- 500
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
MM <- tcrossprod(M)
Y <- diag(1, nrow(MM))


sol <- LOOE(MM,Y)
sol$lambda.min

MMinv <- solveEigen(MM, Y, lambda=sol$lambda.min)
MMinv <- solveEigen(MM, Y, lambda=0.00001)
ans <- MMinv%*%MM
ans <- sol$Ginv%*%MM
diag(ans)



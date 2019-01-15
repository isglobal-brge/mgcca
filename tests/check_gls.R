n <- 500
p <- 100
M <- matrix(rnorm(n*p), nrow=n, ncol=p)
M2 <- cbind(rep(1, n), M)
Y <-  1.2 + 2.4*M[,1] + 1.6*M[,2] - 0.4*M[,5] +
  rnorm(n, 0, 1.2)


ans1 <- glmnet::cv.glmnet(M, Y)
ans1$lambda.min
plot(ans1)

sol <- LOOE(tcrossprod(M),Y)
sol$lambda.min
ans2 <- solveEigen(M2, Y, lambda=sol$lambda.min)

# ans2 <- t(sol$Ginv%*%M)%*%Y  .... this is the same
plot(log(sol$lambdas), sol$looe)

ans3 <- rfunctions::cgls(M2, Y, lambda=sol$lambda.min)

mm <- which(ans1$lambda==ans1$lambda.min)
ans1$glmnet.fit$beta[1:10,mm]
ans2[1:10]
ans3$x[1:10]

plot(ans1)
plot(log(sol$lambdas), sol$looe)


ans4 <- krls(M, Y, whichkernel = "linear")
ans4$lambda
ans4$coeffs[1:5,]

#######################################################
####        OLD
#######################################################

ss <- eigen(MM)

lambda <- 0.06

Lambda <- ss$values
Q <- ss$vectors

LOOE1 <- function(lambda, Lambda, Q, Y){
  W <- 1/(Lambda + lambda)
  Ginv <- rfunctions::crossprodcpp(t(Q), W)
  cte <- Ginv%*%Y
  ans <- sum((cte/diag(Ginv))^2)
  attr(ans, "Ginv") <- Ginv
  ans
}

lambdas <- seq(0.01, 1, length=100)
o <- sapply(lambdas, LOOE, Lambda=Lambda, Q=Q, Y=Y)

# betas

sol <- LOOE(0.006, Lambda, Q, Y)
attr(sol, "Ginv")%*%t(M)%*%Y


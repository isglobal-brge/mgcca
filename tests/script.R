fsim <- function(n=70,p=150,q=150, prop.miss=0.10, ncomp=5, niter=10,
                 tol=1e-5, verbose=FALSE){
  # n = number of individuals
  # p = number of variables in X matrix
  # q = number of variables in Y matrix
  # prop.miss= proportion of missing values
  # ncomp = number of correlation functions
  # niter = maximum number of iterations in mcc algorithm
  # tol = tolerance of mcc algorithm to converge

  ##### simulated data
  u <- matrix(runif(n,-1,1),ncol=1)
  v1 <- matrix(runif(p,-1,1),ncol=1)
  v2 <- matrix(runif(q,-1,1),ncol=1)
  X1 <- u%*%t(v1) + matrix(rnorm(n*p,0,0.1),ncol=p)
  X2 <- u%*%t(v2) + matrix(rnorm(n*q,0,0.1),ncol=q)
  rownames(X1) <- rownames(X2) <- 1:n
  X<-Xori <- list(X1,X2)
  ####### introduce missings at random
  n<- nrow(X[[1]])
  Xmiss <- Xori
  k <- length(Xmiss)
  for (i in 1:k){
    Xi <- X[[i]]
    n <- nrow(Xi)
    p <- ncol(Xi)
    # Xi[sample(1:(n*p),round(prop.miss*n*p))] <- NA
    Xi[sample(1:n,round(prop.miss*n)),] <- NA
    Xmiss[[i]] <- Xi
  }

  #### compare methods (average imputation, mcc and complete rows)

  ## gold
  #st <- system.time(result <- CCA::cc(Xori[[1]], Xori[[2]]))
  result <- mgcca(Xori, nfac = ncomp)
  #Agold <- list(result$xcoef, result$ycoef)
  Agold <- result$A

  ## average imputation
  Ximp <- lapply(Xmiss, function(Xi){
    apply(Xi,2, function(xx){
      m <- mean(xx,na.rm=TRUE)
      ifelse(is.na(xx), m , xx)
    })
  })

#  result <- CCA::cc(Ximp[[1]], Ximp[[2]])
#  Aimp <- list(result$xcoef,result$ycoef)
    result <- mgcca(Ximp, nfac = ncomp)
    Aimp <- result$A


  ## complete cases
  common <- Reduce(intersect,lapply(lapply(Xmiss, na.omit),rownames))
  Xcommon <- lapply(Xmiss, function(xi) xi[common,,drop=FALSE])
  # result <- try(CCA::cc(Xcommon[[1]], Xcommon[[2]]), TRUE)
  # Acommon <- list(result$xcoef,result$ycoef)
  result <- mgcca(Xcommon, nfac = ncomp)
  if (verbose)
    cat("n common=", nrow(Xcommon[[1]]),"\n")
  Acommon <- result$A

  ## mcca
  Xr <- lapply(Xmiss, na.omit)
  result <- mgcca(Xr, nfac=ncomp)
  Ar <- result$A

  ### comparison
  ans <- list()
  for (i in 1:2){
    ans[[i]] <- list(mcc = Agold[[i]] - Ar[[i]],
                     imp = Agold[[i]] - Aimp[[i]],
                     com = Agold[[i]] - Acommon[[i]])
    }
  ans
}

nsim <- 5
out <- list()
for (k in 1:nsim){
  out[[k]] <- fsim(n=100, p=10, q=10, prop.miss=0.1,
             ncomp=10, verbose=TRUE)
}

ans <- matrix(NA, nrow=2, ncol=3)

for (i in 1:2){
  for (j in 1:3){
    temp <- purrr::map(purrr::map(out,i),j)
    ans[i,j] <- sum(Reduce("+", lapply(temp, function(x) x^2)))
  }
}
colnames(ans) <- c("mcc", "impute", "complete")
rownames(ans) <- c("X", "Y")
ans



#
#  Un set de datos
#

set.seed(12345)
n <- 500
p <- 140
q <- 75
prop.miss <- 0.1
ncomp <- 2

u <- matrix(runif(n,-1,1),ncol=1)
v1 <- matrix(runif(p,-1,1),ncol=1)
v2 <- matrix(runif(q,-1,1),ncol=1)
X1 <- u%*%t(v1) + matrix(rnorm(n*p,0,0.1),ncol=p)
X2 <- u%*%t(v2) + matrix(rnorm(n*q,0,0.1),ncol=q)
rownames(X1) <- rownames(X2) <- 1:n
colnames(X1) <- paste0("var", 1:ncol(X1))
colnames(X2) <- paste0("var", 1:ncol(X2))
X<-Xori <- list(X1,X2)

####### introduce missings at random
n<- nrow(X[[1]])
Xmiss <- Xori
k <- length(Xmiss)
for (i in 1:k){
  Xi <- X[[i]]
  n <- nrow(Xi)
  p <- ncol(Xi)
  # Xi[sample(1:(n*p),round(prop.miss*n*p))] <- NA
  Xi[sample(1:n,round(prop.miss*n)),] <- NA
  Xmiss[[i]] <- Xi
}

###### Get tables with missing rows
Xr <- lapply(Xmiss, na.omit)


###### Libraries
library(mgcca)
library(RGCCA)
library(rfunctions)
library(CCA)
library(omicade4)

system.time(result1 <- mgcca(Xori, nfac = ncomp))
system.time(result2 <- mgcca(Xori, nfac = ncomp, method="penalized",
                             lambda=c(0.5, 0.1)))
system.time(result3 <- mgcca(Xori, nfac = ncomp, method="geninv"))



C <- matrix(c(0,0,1,
              0,0,1,
              1,1,0), 3, 3)

Xg <- list(Xori[[1]], Xori[[2]], cbind(Xori[[1]], Xori[[2]]))
system.time(result4 <- rgcca(Xg, C=C, tau = rep(0,3),
                 scheme = "factorial", ncomp=c(2,2,2), verbose=FALSE))

system.time(result5 <- mcia(list(t(Xori[[1]]), t(Xori[[2]]))))


microbenchmark(
  mgcca = res1 <- mgcca(Xori, nfac = ncomp),
  gcca = res2 <- gcca(Xori, nfac = ncomp),
  times = 100
)


# some checks
a <- matrix(rnorm(10000), 100, 100)
aa <- crossprod(a)
ee <- eigen(aa)
Q <- ee$vectors
L <- ee$values
lambda <- 0.2

ans1 <- solve_penal(Q, L, lambda)
ans2 <- chol2inv(chol(aa + diag(ncol(aa))*lambda))
ans3 <- .C("solve_penal", Q, L, lambda, as.integer(nrow(Q)))

system.time(

ff <- function(x){
  a1<-RConics::adjoint(x)
  a2 <- det(x)
  a3 <- a1/a2
  a3
}
system.time(o2 <- ff(aa))
system.time(o2 <- solve(aa))


ans1[1:5, 1:5]
ans2[1:5, 1:5]

LOOE <- function(lambda, Q, L){
  num <- getC(Q, L, lambda)
  dem <- getGinv(Q, L, lambda)
  dd <- num/dem
  ans <- sum(diag(dd)^2)
  ans
}

LOOE <- function(lambda, Q, L){
  ans <- getC(Q, L, lambda)
  sum(diag(ans))
}



getC <- function(Q, L, lambda){
  inv <- 1/(L+lambda)
  QI <- mult_Xw(Q, inv)
  QIQ <- tcrossprod(QI, Q)
  QIQ
}

getGinv <- function(Q, L, lambda){
  ans <- sapply(1:length(L), getGinv.i, Q=Q, L=L, lambda=lambda)
  ans
}

getGinv.i <- function(j, Q, L, lambda){
  ans <- sapply(1:length(L), getGinv.j, Q=Q, L=L, lambda=lambda, j=j)
  ans
}

getGinv.j <- function(i, j, Q, L, lambda){
  v <- (Q[i,]*Q[j,])/(L+lambda)
  ans <- sum(v)
  ans
}

system.time(ans3 <- outer(1:length(L), 1:length(L), getGinv.j,
                          Q=Q, L=L, lambda=0.2))

solve_penal <- function(Q, L, lambda){
  # Q square matrix (nxn)
  # L vector (n)
  # lambda ()
  n <- nrow(Q)
  ans <- matrix(NA, n, n)
  for (i in 1:n){
    for (j in i:n) {                       # nota que va desde i porque es simetrica
      ans.ij <- (Q[i,]*Q[j,])/(L+lambda)   # nota que el bucle para hacer
      ans[i,j] <- ans[j, i] <- sum(ans.ij) # la suma en L no la pongo porque
    }                                      # R es vectorial
  }
  ans
}



# check new algorithm

ll <- estim.regul(X[[1]], X[[2]])

mod <- mgcca(X, method="penalized", lambda=c(1,0.1))
i <- 1
xk <- mult_Xw(t(X[[i]]), diag(K[[i]]))
M <- xk%*%X[[i]]

ee <- eigen(M)
Q <- ee$vectors
L <- ee$values
ss <- seq(0.05, 1, 0.05)
ans <- sapply(ss, LOOE, Q=Q, L=L)
ss[which.min(ans)]

# plot

plot(result4$Y[[3]], type="n")
text(result4$Y[[3]], rownames(Xori[[1]]))

plot(result1$Y, type="n")
text(result1$Y, rownames(Xori[[1]]))

plot(result1$corsY[[1]], type="n")
text(result1$corsY[[1]], rownames(Xori[[1]]))

# Check diag(K)%*%M%*%diag(K)
A <- matrix(rnorm(25), 5, 5)
w <- round(seq(1,2, length=5),1)
res1 <- diag(w)%*%A%*%diag(w)
res2 <- mult_wXw(A,w)
identical(res1,res2)

res3 <- crossprod(A, diag(w))
res4 <- mult_Xw(t(A),w)
identical(res3,res4)


# check with CCA

mod1 <- mgcca(X, scores=TRUE)
mod2 <- cc(X1, X2)
plotInds(mod1, print.labels = TRUE)
plt.indiv(mod2, 1, 2)

data(nutrimouse, package="CCA")
X <- as.matrix(nutrimouse$gene)
Y <- as.matrix(nutrimouse$lipid)
group <- nutrimouse$genotype
X.s <- scale(X)
Y.s <- scale(Y)
lambda.est <- estim.regul(X.s, Y.s, plt=FALSE)
lambda <- c(lambda.est[1], lambda.est[2])
res.rcc <- rcc(X.s, Y.s, lambda[1], lambda[2])
res.mgcca <- mgcca(list(X.s, Y.s), scale=FALSE, method="penalized",
                   lambda=c(0.13,0.13))

plotInds(res.mgcca, print.labels = TRUE)
plotInds(res.mgcca, group=group, print.labels = TRUE)
plt.indiv(res.rcc, 1, 2)


# SVD big data

a <- matrix(rnorm(10000), 100, 100)
aa <- crossprod(a)
system.time(corpcor::fast.svd(a))
system.time(RSpectra::svds(a, k=5))
system.time(svd(a))

system.time(aa <- solve(a))
system.time(aaa <- rfunctions::cgls(a, diag(ncol(a)))$x)
system.time(aaaa <- rfunctions::solveEigen(a, diag(ncol(a)))$x)

# check inverse

a <- matrix(rnorm(50000000), 5000, 10000)
a <- matrix(rnorm(50000000), 5000, 10000)
aa <- crossprod(a)
aaa <- aa + 0.2*diag(ncol(aa))

system.time(ans1 <- chol2inv(chol(aaa)))
system.time(ans2 <- bigLmFit(diag(nrow(aaa)), aaa, chunk=100, byrow=TRUE))

a1 <- inv2(aa, lambda=0.000001)
sum((a%*%a1%*%t(a) - diag(nrow(a)))**2)



inv2 <- function(x, lambda){
  y <- diag(nrow(x))*lambda
  M <- x + y
  ans <- chol2inv(chol(M))
  ans
}

system.time(ans1 <- inv2(aa, 0.3))
ee <- eigen(aa)
Q <- ee$vectors
L <- ee$values
system.time(ans2 <- getC(Q, L, 0.3))


getC <- function(Q, L, lambda){
  Xw <- mult_Xw(Q, 1/(L+lambda))
  ans <- tcrossprod(Xw, Q)
  ans
}

getC1 <- function(x, lambda){
  eig <- eigen(x, symmetric=TRUE)
  Q <- eig$vectors
  L <- eig$values
  Xw <- mult_Xw(Q, 1/(L+lambda))
  ans <- tcrossprod(Xw, Q)
  ans
}


LOOE <- function(lambda, x) {
  inv <- inv2(x, lambda)
  looe <- inv - diag(nrow(inv))
  ans <- sum(looe^2)
  ans
}

LOOE <- function(lambda, Q, L) {
  inv <- getC(Q, L, lambda)
  looe <- inv - diag(nrow(inv))
  ans <- sum(looe^2)
  ans
}

ee <- eigen(aa)
Q <- ee$vectors
L <- ee$values

M <- getC(Q, L, 0.2)
(MM%*%M)[1:5,1:5]

ss <- seq(0.05, 4, 0.01)
ans <- sapply(ss, LOOE, x=aa)
ans

ok <- ss[which.min(ans)]

(aa%*%getC(Q,L, 0.1))[1:5,1:5]



############################

M <- cgls(MM, diag(ncol(MM)), lambda=0.4)$x
Minv <- MM%*%M
Minv[1:5,1:5]


ee <- eigen(M, symmetric = TRUE)
Q <- ee$vectors
L <- ee$values
ss <- seq(0.05, 1, 0.05)
ans <- sapply(ss, LOOE, Q=Q, L=L)
ss[which.min(ans)]



b1 <- solve(aa)
b2 <- cgls(aa, y, lambda=0)$x
b3 <- RIDGEsigma(aa, lam=c(0.001,0.2))

max(b4-b1)
max(b2-b1)
max(b3-b1)


b1 <- solve(M)
b2 <- speedglm::speedlm.fit(c(1,nrow(M)-1), M)

# check ginverse

a <- matrix(rnorm(1000000), 100, 100)
a <- crossprod(a)
matrixcalc::is.positive.definite(a)

system.time(b1 <- ginv(a))
system.time(b2 <- rfunctions::geninv(a))

max(b1-b2)



#
# NUTRI
#

data(nutrimouse)
X <- as.matrix(nutrimouse$gene[,1:5])
Y <- as.matrix(nutrimouse$lipid)

X <- scale(X)
Y <- scale(Y)

C <- matrix(c(0,0,1,
              0,0,1,
              1,1,0), 3, 3)

Xg <- list(X, Y, cbind(X,Y))

res.cc <- cc(X,Y)
res.mgcca <- mgcca(list(X,Y), scale=FALSE, scores=TRUE)
res.rgcca <- rgcca(Xg, C=C, tau = rep(0,3),
                   ncomp=rep(2,3), verbose=FALSE,
                   scheme="factorial", scale=FALSE)


max(abs(res.mgcca$scores[[1]]) - abs(res.cc$scores$xscores[,1:2]))
max(abs(res.mgcca$scores[[2]]) - abs(res.cc$scores$yscores[,1:2]))

max(abs(res.mgcca$Y) - abs(res.rgcca$Y[[3]]))


#
# Russet
#
library(RGCCA)
data(Russett)
X_agric <- as.matrix(Russett[,c("gini","farm","rent")])
X_ind <- as.matrix(Russett[,c("gnpr","labo")])
X_polit <- as.matrix(Russett[ , c("inst", "ecks", "death",
                                  "demostab", "dictator")])
A = list(X_agric, X_ind, X_polit, cbind(X_agric, X_ind, X_polit))

#Define the design matrix (output = C)
C <- matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0), 4, 4)
result.rgcca <- rgcca(A, C, tau = c(0, 0, 0, 0), ncomp = rep(2, 4),
                      scheme = "factorial", scale = TRUE, verbose=FALSE)
lab <- as.vector(apply(Russett[, 9:11], 1, which.max))
names(lab) <- rownames(Russett)

plot(result.rgcca$Y[[4]][, 1], result.rgcca$Y[[4]][, 2], col = "white",
     xlab = "Global Component 1", ylab = "Global Component 2")
text(result.rgcca$Y[[4]][, 1], result.rgcca$Y[[4]][, 2], rownames(Russett),
     col = lab, cex = .7)
grid()
abline(h=0)
abline(v=0)

library(microbenchmark)
microbenchmark(
  mgcca = res.mgcca <- mgcca(list(X_agric, X_ind, X_polit), scores=FALSE),
  gcca = res.gcca <- gcca(list(X_agric, X_ind, X_polit), scores=FALSE)
)


plotVar(res.mgcca, cor.thres = 0, legend.show = FALSE, tcolors = c("red", "blue", "darkgreen"),
        lab.cex = 3, point.cex = 2)
plotInd(res.mgcca, group=as.factor(lab))


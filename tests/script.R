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
n <- 1000
p <- 200
q <- 500
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
                             lambda=c(0.1, 0.1)))

system.time(result3 <- mcia(list(t(Xori[[1]]), t(Xori[[2]]))))



C <- matrix(c(0,0,1,
              0,0,1,
              1,1,0), 3, 3)

Xg <- list(Xori[[1]], Xori[[2]], cbind(Xori[[1]], Xori[[2]]))
system.time(result4 <- rgcca(Xg, C=C, tau = rep(0,3),
                 scheme = "factorial", ncomp=c(2,2,2), verbose=FALSE))


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

data(nutrimouse)
X <- as.matrix(nutrimouse$gene)
Y <- as.matrix(nutrimouse$lipid)
group <- nutrimouse$genotype
X.s <- scale(X)
Y.s <- scale(Y)
lambda.est <- estim.regul(X.s, Y.s, plt=FALSE)
lambda <- c(lambda.est[1], lambda.est[2])
res.rcc <- rcc(X.s, Y.s, lambda[1], lambda[2])
res.mgcca <- mgcca(list(X.s, Y.s), scale=FALSE, method="penalized",
                   lambda=lambda)

plotInds(res.mgcca, print.labels = TRUE)
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

a <- matrix(rnorm(10000), 100, 100)
aa <- crossprod(a)
y <- diag(nrow(a))

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


res.mgcca <- mgcca(list(X_agric, X_ind, X_polit), scores=TRUE)


plotVar(res.mgcca, cor.thres = 0, legend.show = FALSE, tcolors = c("red", "blue", "darkgreen"),
        lab.cex = 3, point.cex = 2)
plotInd(res.mgcca, group=as.factor(lab))


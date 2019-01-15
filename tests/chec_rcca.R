library(CCA)
library(mgcca)
data(nutrimouse, package = "CCA")
X1 <- as.matrix(nutrimouse$gene)
X2 <- as.matrix(nutrimouse$lipid)
group <- nutrimouse$genotype

ll <- cv_gcca(X1, X2)
mm <- mgcca(list(X1, X2), method="penalized",
            lambda=c(0.75,1))
plotInds(mm, group=group)

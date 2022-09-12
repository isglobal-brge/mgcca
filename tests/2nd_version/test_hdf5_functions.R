
library(mgcca)
library(BigDataStatMeth)
library(rhdf5)

# Working directory
setwd("/Users/mailos/DOCTORAT_Local/mgcca")

load("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca/bd/multiassayexperiment.rda")
load("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca/bd/multiassayexperiment.imputed.rda")

## 3) Obtain a list of tables (matrices), one for each assay
tables_list <- getTables(multiassayexperiment.imputed)

# The sample ID is the first 12 characters from the sample name
for (assay in 1:length(tables_list)) {
  rownames(tables_list[[assay]]) <- substr(rownames(tables_list[[assay]]),1,12)
}

# Subset data
table_list.subset <- lapply(tables_list, function(x) x[,1:1000])




## Test hdf5 functions


#..## Proves prèvies a execució mgcca.... però només son probes, la part bona es a parrir de mgcca_bd

tables_list.subset <- getTables(multiassayexperiment.imputed)

# Get data in datasets inside an hdf5 data file
tables_list.subset <- mgcca::getTables_hdf5(multiassayexperiment.imputed, "tmp/gettables.hdf5", overwriteFile = T, overwriteDataset = T)
bdNormalize_hdf5("tmp/gettables.hdf5", "MGCCA_IN", "ACC_Methylation-20160128", wsize = 100)

h5fsvd = H5Fopen("tmp/gettables.hdf5")
matnorm <- h5fsvd$NORMALIZED$MGCCA_IN$`ACC_Methylation-20160128`
colnames(matnorm) <- t(h5fsvd$MGCCA_IN$`.ACC_Methylation-20160128_dimnames`$`2`)
rownames(matnorm) <- t(h5fsvd$MGCCA_IN$`.ACC_Methylation-20160128_dimnames`$`1`)
h5closeAll()

dscalat <- scale(multiassayexperiment.imputed[[2]])
rownames(dscalat) <- NULL
colnames(dscalat) <- NULL

dscalat[1:5,1:5]
matnorm[1:5,1:5]
all.equal(dscalat[1:5,1:5],matnorm[1:5,1:5], check.attributes = FALSE)
dim(dscalat)
dim(matnorm)




# devtools::reload(pkgload::inst("BigDataStatMeth"))
# devtools::reload(pkgload::inst("mgcca"))
#

# Code form mgcca_hdf5


library(mgcca)
library(BigDataStatMeth)
library(rhdf5)

# Working directory
setwd("/Users/mailos/DOCTORAT_Local/mgcca")

#. Dades sense imputar - no cal carregar-les .# load("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca/bd/multiassayexperiment.rda")
load("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca/bd/multiassayexperiment.imputed.rda")

# devtools::reload(pkgload::inst("BigDataStatMeth"))
# devtools::reload(pkgload::inst("mgcca"))

tables_list <- getTables(multiassayexperiment.imputed)
# Modifiquem els noms dels ids i els passem a només 12 caràcters que en realitat es el id de la mostra, no els
for ( i in 1:length(tables_list)) {
  rownames(tables_list[[i]]) <- substr(rownames(tables_list[[i]]), 1, 12)
}

doubleExp <- list("ACC_RNASeq2GeneNorm-20160128"  = t(tables_list[[1]]), "ACC_Methylation-20160128" = t(tables_list[[2]]))
multiassayexperiment.imputed <- MultiAssayExperiment(experiments=doubleExp, colData = colData(multiassayexperiment.imputed))



# REDUIM EL TAMANY DE LES DADES PER PODER FER PROVES COM CAL....
tables_list.subset <- mgcca::getTables_hdf5(multiassayexperiment.imputed[1:100, ,], "tmp/gettables.hdf5", overwriteFile = T, overwriteDataset = T)
# Totes le dades....
#..# tables_list.subset <- mgcca::getTables_hdf5(multiassayexperiment.imputed, "tmp/gettables.hdf5", overwriteFile = T, overwriteDataset = T)

devtools::reload(pkgload::inst("BigDataStatMeth"))
devtools::reload(pkgload::inst("mgcca"))

hh <- mgcca_hdf5(x = "tmp/gettables.hdf5",
                 filename = "tmp/gettables.hdf5",
                 group = "MGCCA_IN",
                 datasets = c("ACC_Methylation-20160128", "ACC_RNASeq2GeneNorm-20160128"),
                 method="penalized", lambda=c(0.75,1))

hh <- mgcca_hdf5(x = "tmp/gettables.hdf5",
                 filename = "tmp/gettables.hdf5",
                 group = "MGCCA_IN",
                 datasets = c("ACC_Methylation-20160128", "ACC_RNASeq2GeneNorm-20160128"),
                 method="penalized", lambda=c(1, 0.75), scores = T)




mat1.data <- c(1,2,3,4,5,6,7,8,9)
mat1 <- matrix(mat1.data,nrow=3,ncol=3,byrow=TRUE)


BigDataStatMeth::bdCreate_hdf5_matrix_file("pepet.hdf5", object = mat1, group = "aa", dataset = "aaa")
BigDataStatMeth::bdCreateGroup_hdf5("pepet.hdf5", "tmp/test")



# Code form mgcca and mgcca_bd

a <- mgcca(table_list.subset, method="penalized", lambda=c(0.75,1))
b <- mgcca_bd(table_list.subset, method="penalized", lambda=c(0.75,1))


# Plot data
group <- factor(multiassayexperiment$vital_status)

par(mfrow = c(2,1))
  plotInds(a, group=group, pos.leg = "bottomleft")
  plotInds(b, group=group, pos.leg = "bottomleft")
par(mfrow = c(1,1))





library(microbenchmark)
res <- microbenchmark( mgcca(table_list.subset, method="penalized", lambda=c(0.75,1)),
                       mgcca_bd(table_list.subset, method="penalized", lambda=c(0.75,1)),
                       times = 1)
res

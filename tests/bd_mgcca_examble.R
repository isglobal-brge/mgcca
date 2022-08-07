
library(mgcca)
library(R.utils)

# libraries needed by TCGA data
library(curatedTCGAData)
library(TCGAutils)

# Libraries needed by mgcca package
library(MultiAssayExperiment)
library(RSpectra)
library(impute)


## 1) Download TCGA data

# Available cancer codes and assays
curatedTCGAData(diseaseCode = "*", assays = "*", dry.run = TRUE)
# MultiAssayExperiment created with two assays: RNA-seq and Methylation
#multiassayexperiment <- curatedTCGAData("ACC", c("RNASeq2GeneNorm","Methylation", "miRNASeqGene","RPPAArray" ), FALSE)
multiassayexperiment <- curatedTCGAData("ACC", c("RNASeq2GeneNorm","Methylation" ), FALSE)
multiassayexperiment


# Now let's check how this MultiAssayExperiment looks like.


class(multiassayexperiment)
multiassayexperiment
# Methylation data
assay(multiassayexperiment[[1]])[1:5,1:5]
dim(assay(multiassayexperiment[[1]]))
# RNA-seq data
assay(multiassayexperiment[[2]])[1:5,1:5]
dim(assay(multiassayexperiment[[2]]))
# We can see there are some missing values (NA and/or 0)
table(is.na(assay(multiassayexperiment[[1]])))
table(is.na(assay(multiassayexperiment[[2]])))

## 2) Imputation

# Arguments:
# subset: MultiAssayExperiment
# method: knn for K-Nearest Neighbors or
# k: Number of neighbors to be used in the imputation
# rowmax: maximum percent missing data allowed in any row (default 50%). For any rows with more than rowmax% missing are imputed using the overall mean per sample.
# colmax: maximum percent missing data allowed in any column (default 80%). If any column has more than colmax% missing data, the program halts and reports an error.
# remove.col: To remove columns with more than colmax% missing data.
# impute.zero: To impute 0 values.

multiassayexperiment.imputed <- impute(multiassayexperiment, method = "knn", k = 10,
                                       rowmax = 0.5, colmax = 0.8, remove.col = TRUE, impute.zero = TRUE)

## 3) Obtain a list of tables (matrices), one for each assay


tables_list <- getTables(multiassayexperiment.imputed)

# We assess the class, length, data and dimensions of the list of tables.
class(tables_list)
length(tables_list)
tables_list[[1]][1:5,1:5]
tables_list[[2]][1:5,1:5]
dim(tables_list[[1]])


######
######
######
####   ESTIC PER AQUÍ FENT LES MEVES PROVETES....
####   Que faig amb els fitxers??? ho executo sobre fitxer ????
####   si ho faig sobre fitxer --> Escric cada una de les matrius al fitxer i enlloc de treballar sobre
####   llistes treballo sobre datasets???
####   Cada element de la llista es un dataset --> ho podria guardar d'alguna forma que jo sàpiga quin es quin....
####   potser podria guardar amb la estructura :
####      multiassayexperiment <nom del multiassay>:
####         Dataset[1] :<<atribut : nom de la llista>>
####                     Dataset (matriu de dades)
####         Dataset[2] :<<atribut : nom de la llista>>
####                     Dataset (matriu de dades)
####         Dataset[3] :<<atribut : nom de la llista>>
####                     Dataset (matriu de dades)
####         ......


####         ##. BDSM .## Convert MultiassayExperimetn to HDF5
####         library(HDF5Array)
####
####         # We use rhdf5 package to write data to disk, we also can do dat with BDSM package
####         testh5 <- tempfile(fileext = ".h5")
####         writeHDF5Array(smallMatrix, filepath = testh5, name = "smallMatrix",
####                        with.dimnames = TRUE)
####         ##. BDSM .##


######
######
######



## 4) Additional pre-processing (depending on your data)
# Before running `mgcca`, we need to make some changes that depends on the type of data you have.

# Sample names


# The sample ID is the first 12 characters from the sample name
for (assay in 1:length(tables_list)) {
  rownames(tables_list[[assay]]) <- substr(rownames(tables_list[[assay]]),1,12)
}


# Tables containing character numbers

# In our case, methylation data seems to be of character, instead of numeric, and although the imputation can be done with that type of data this will give an error in `mgcca`, because it expects numeric matrices. Thus, let's change this matrix to numeric values.


# Check if it is character
tables_list[[1]][1,1]
class(tables_list[[1]][1,1])

# Use our in-house function
tables_list <- matrix.chr2num(tables_list, index = 1)

# Confirm it's numeric now
tables_list[[1]][1,1]
class(tables_list[[1]][1,1])

# 4) MGCCA


class(tables_list[[1]])
table_list.subset <- lapply(tables_list, function(x) x[,1:100])
lapply(table_list.subset, function(x) x[1:5,1:5])

# Run mgcca, select the grouping variable for the plot, and finally plot the result.

# mgcca
mgcca.result <- mgcca(table_list.subset, method="penalized", lambda=c(0.75,1))


######    BDSM    #######

# mgcca BDSM
mgcca.result.bd <- mgcca_bd(table_list.subset, method="penalized", lambda=c(0.75,1))

######    BDSM    #######


mtime <- microbenchmark::microbenchmark(mgcca.result <- mgcca(table_list.subset, method="penalized", lambda=c(0.75,1)),
                                        mgcca.result.bd <- mgcca_bd(table_list.subset, method="penalized", lambda=c(0.75,1)),
                                        times = 5)

print(summary(mtime)[, c(1:7)],digits=3)


# grouping variable from the original MultiAssayExperiment
# In this case, as example, the vital status of the patients with ACC cancer.
group <- factor(multiassayexperiment$vital_status)
# miRNA cluster
#group <- factor(multiassayexperiment$miRNA.cluster)

#plot
plotInds(mgcca.result, group=group, pos.leg = "topright")

######    BDSM    #######
#plot
plotInds(mgcca.result.bd, group=group, pos.leg = "topright")
######    BDSM    #######


# 5) Mgcca for Genomic Ranges

# You can also make a `mgcca` for a specific genomic coordinate that you are interested in. To do so, we will first subset the features from each assay (RNA-seq and methylation) that belong to that particular genomic coordinate.
# AS an example, we will use the coordinates for two genes:
# IGF2 --> a prototypical tumor-progenitor gene, which is the most commonly over-expressed gene in human adrenocortical carcinoma.
# MEN1 --> A tumor suppressor gene located at 11q13, also related to this cancer.
# TP53 --> The most common tumor supressor gene in all cancers.
# We will use both gene coordinates + 2kb up and downstream (in order to obtain CpGs next to it).
#
# NOTE 1: we need at least 2 genes, because (I think) imputation cannot be done if there is only 1 row (gene).
# NOTE 2: The first thing to do is to find the genes names (and cpgs) belonging to that particular genomic coordinates, but in this case we already know their names. That means that genomic coordinates can be given without knowing the actual genes, because if you know the actual genes you just create a list with its names (see below).

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

### 1) Save the genomic coordinates that we are interested in
# MEN1 --> chr11: 64570986-64578766 +-2kb --> 64568986-64580766

chr.coords <- GRanges("11:64550986-64578766")

### 2) Annotate each assay independently to obtain the genomic coordinates for all features, in case it is not already done

### 2.1) Methylation annotation
# We already have the genomic coordinates for each CpGs in the rowData dataframe
rowData(multiassayexperiment[[1]])

### 2.2) RNA-seq annotation
# It will split the assay in two: ranged RangedSummarizedExperiment containing the features with existing genomic coordinates, and unranged SummarizedExperiment, not containing them
mae.annotated <- TCGAutils::symbolsToRanges(multiassayexperiment)

### 3) Remove the non-used assays (unranged SummarizedExperiment)

mae.annotated[[3]] <- NULL

### 4) Obtain a subset of features (CpGs and genes) that fall within the interested genomic coordinates

mae.subset <- genCoords(multiassayexperiment = mae.annotated , gen.coords = chr.coords,
          meth.index = 1, col.name = "Genomic_Coordinate")

# mae.subset <- genCoords(multiassayexperiment = mae.annotated , gen.coords = chr.coords,
#           type.omics = c("methylation","rnaseq"), col.name = "Genomic_Coordinate")


# Now we can make the `mgcca` with the MultiAssayExperiment subset by genomic coordinates


### 1) Imputation

mae.subset.imputed <- impute(mae.subset, method = "knn", k = 10,
                                       rowmax = 0.5, colmax = 0.8, remove.col = TRUE, impute.zero = TRUE)

### 2) List of Matrices

tables_list.subset <- getTables(mae.subset.imputed)

### 3) Reassign sample names

for (assay in 1:length(tables_list.subset)) {
  rownames(tables_list.subset[[assay]]) <- substr(rownames(tables_list.subset[[assay]]),1,12)
}

### 4) Convert Character matrix (methylation) to Numeric matrix

tables_list.subset <- matrix.chr2num(tables_list.subset, index = 1)

### 5) mgcca

# mgcca
mgcca.result <- mgcca(tables_list.subset, method="penalized", lambda=c(0.75,1))


######    BDSM    #######

## ................... INICI PROVES ................... ##

#..## Proves prèvies a execució mgcca.... però només son probes, la part bona es a parrir de mgcca_bd

tables_list.subset <- getTables(mae.subset.imputed)

# Get data in datasets inside an hdf5 data file
tables_list.subset <- mgcca::getTables_hdf5(mae.subset.imputed, "tmp/gettables.hdf5")
Normalize_hdf5("tmp/gettables.hdf5", "MGCCA_IN", "ACC_Methylation-20160128", wsize = 3)

h5fsvd = H5Fopen("tmp/gettables.hdf5")
matnorm <- h5fsvd$NORMALIZED$MGCCA_IN$`ACC_Methylation-20160128`
h5closeAll()

dscalat <- scale(tables_list.subset[[1]])
rownames(dscalat) <- NULL
colnames(dscalat) <- NULL

dscalat[1:5,1:5]
matnorm[1:5,1:5]
all.equal(dscalat[1:5,1:5],matnorm[1:5,1:5], check.attributes = FALSE)
dim(dscalat)
dim(matnorm)

## ................... FI PROVES ................... ##



# mgcca
mgcca.result <-  mgcca_bd(tables_list.subset, method="penalized", lambda=c(0.75,1))

######    BDSM    #######

mtime <- microbenchmark::microbenchmark(mgcca.result <- mgcca(tables_list.subset, method="penalized", lambda=c(0.75,1)),
                                        mgcca.result_bd <- mgcca_bd(tables_list.subset, method="penalized", lambda=c(0.75,1)),
                                        times = 15)

print(summary(mtime)[, c(1:7)],digits=3)





# grouping variable from the original MultiAssayExperiment
# In this case, as example, the vital status of the patients with ACC cancer.
group <- factor(multiassayexperiment$vital_status)
# miRNA cluster
#group <- factor(multiassayexperiment$miRNA.cluster)

#plot
plotInds(mgcca.result, group=group, pos.leg = "bottomleft")


######    BDSM    #######

# plot
plotInds(mgcca.result_bd, group=group, pos.leg = "bottomleft")

######    BDSM    #######

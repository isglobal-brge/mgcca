####
####  BENCHMARK mgcca
####
####  Date : 25/01/2022
####  Author : Dolors Pelegrí
####
####  Description : Benchmark for mgcca using :
####      - R-functions with
####          . R-Objects
####      - BigDataStatMeth library with :
####          . R-Objects
####          . DelayedArray objects
####          . HDF5 data files
####

library(microbenchmark)

library(BigDataStatMeth)
library(DelayedArray)
library(mgcca)

library(ggthemes)
library(gtable)
library(gridExtra)
library(ggplot2)
library(ggpubr)


library(mgcca)
library(R.utils)

# libraries needed by TCGA data
library(curatedTCGAData)
library(TCGAutils)

# Libraries needed by mgcca package
library(MultiAssayExperiment)
library(RSpectra)
library(impute)

setwd("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca") #..LOCAL..#
tipusdades <- c('R ','DelayedArray ') # tipus execució paral·lel, sequencial



## 1) Download TCGA data
# multiassayexperiment <- curatedTCGAData("ACC", c("RNASeq2GeneNorm","Methylation" ), dry.run = FALSE, version = "2.0.1")
# save(multiassayexperiment, file = "bd/multiassayexperiment.rda")
load("bd/multiassayexperiment.rda")

## 2) Impute TCGA data
# multiassayexperiment.imputed <- impute(multiassayexperiment, method = "knn", k = 10,
#                                        rowmax = 0.5, colmax = 0.8, remove.col = TRUE, impute.zero = TRUE)
# save(multiassayexperiment.imputed, file = "bd/multiassayexperiment.imputed.rda")
load("bd/multiassayexperiment.imputed.rda")

## 3) Obtain a list of tables (matrices), one for each assay
tables_list <- getTables(multiassayexperiment.imputed)

# The sample ID is the first 12 characters from the sample name
for (assay in 1:length(tables_list)) {
  rownames(tables_list[[assay]]) <- substr(rownames(tables_list[[assay]]),1,12)
}


repet <- 1
ncores <- 3

resultsmgca.df <- data.frame( expr = character(),min=numeric(),lq=numeric(),mean=numeric(),median=numeric(),uq=numeric(),
                              max=numeric(),tipusop=character(),tipusdades=character(),bloc_size=numeric(),M=numeric(),
                              K=numeric(), N=numeric(), ncores = numeric(), nexec=numeric())

for ( i in seq(500,2000, by = 500) )
{

  # Subset data
  table_list.subset <- lapply(tables_list, function(x) x[,1:i])

  res <- microbenchmark( mgcca(table_list.subset, method="penalized", lambda=c(0.75,1)),
                         mgcca_bd(table_list.subset, method="penalized", lambda=c(0.75,1)),
                         times = repet, unit = "s")

  resdata <- as.data.frame(summary(res)[, c(1:7)])

  resdata <- cbind(resdata, rbind('R-code','BigDataStatMeth'),tipusdades = 'memory', bloc_size = i,
                   M=dim(table_list.subset[[1]])[1], K=dim(table_list.subset[[1]])[2], N=0, ncores, repet)

  colnames(resdata)[8] <- c('tipusop')
  resultsmgca.df <- rbind(resultsmgca.df,resdata)

}


filename <- "MGCCA-test"
write.csv(resultsmgca.df, paste0("results/",filename,".csv"))

#..#  Read data from file (if needed)
#..# resultsmgca.df <- read.csv(paste0("results/",filename,".csv"))

# Grafiquem els resultats
# lcols <- metagen::cbbPalette[c(2,3)]
linetype = rep(c('solid', 'dashed'))

fname <- paste0("plots/", filename, ".png")
png(fname, width = 10, height = 8, units = 'in', res = 900, pointsize = 4)
p2 <- ggplot(resultsmgca.df, aes( x = bloc_size,
                                  y = mean )) +
  geom_line(aes(color= tipusop ) ) +
  labs(# title = "MGCCA - implementation R and BigDataStatMeth",
    x = "Dimension",
    y = "Time (s)") +
  labs(color=c('Function','Data Type'))+
  theme(text = element_text( size = 14),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title=element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal"
  ) +
  geom_rangeframe() +
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  scale_linetype_manual(name = 'Data Type', values = linetype)
print(p2)

dev.off()




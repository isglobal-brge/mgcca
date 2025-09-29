
library(mgcca)

load("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca/bd/multiassayexperiment.rda")
load("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/Benchmarks/mgcca/bd/multiassayexperiment.imputed.rda")

## 3) Obtain a list of tables (matrices), one for each assay
tables_list <- getTables(multiassayexperiment.imputed)

# The sample ID is the first 12 characters from the sample name
for (assay in 1:length(tables_list)) {
  rownames(tables_list[[assay]]) <- substr(rownames(tables_list[[assay]]),1,12)
}

# Subset data
table_list.subset <- lapply(tables_list, function(x) x[,1:100])


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

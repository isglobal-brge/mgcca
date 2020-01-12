#' Subset a MultiAssayExperiment object by specific genomic coordinates
#' @param gen.coords
#' @param meth.index
#' @param col.name

# gen.coords must be a GRanges object
# meth.index is the index where the methylation data (or a SummarizedExperiment non Ranged) is located inside the mae.

genCoords <- function(multiassayexperiment, gen.coords, meth.index, col.name){
  
  # Check that the input is a MultiAssayExperiment
  if (!class(multiassayexperiment) == "MultiAssayExperiment")
    stop("Input must be a 'MultiAssayExperiment' object \n")
  
  # Check that gen.coords is a GRanges object
  if (!class(gen.coords) == "GRanges")
    stop("'gen.coords' must be a 'GRanges' object \n")

  subset.names = list()
  
  for (assay in 1:length(multiassayexperiment)) {
    
    # Check which type of omics data is each assay
    # No RangedSummarizedExperiment (e.g. Methylation)
    if (assay == meth.index) {
      # Obtain CpGs names
      subset.names[[names(experiments(multiassayexperiment)[assay])]] <- norangedCoords(multiassayexperiment[[assay]], gen.coords, col.name = "Genomic_Coordinate")
    }
    # RangedSummarizedExperiment (e.g. RNA-seq, miRNA, proteomics, etc)
    else {
      # Obtain feature names
      subset.names[[names(experiments(multiassayexperiment)[assay])]] <- rangedCoords(multiassayexperiment[[assay]], gen.coords)      
    }
  }

  # MultiAssayExperiment subset using feature names obtained above with our genomic coordinate
  mae.subset <- subsetByRow(multiassayexperiment, subset.names)
  mae.subset
  
}

# Subset of features within the specific genomic coordinates for a RangedSummarizedExperiment object

rangedCoords <- function(rangedsummarizedexperiment, gen.coords){
  
  # Check if there is rowData (annotation)
  if (length(rowRanges(rangedsummarizedexperiment)) == 0)
    stop("There is no rowRanges data in the given `RangedSummarizedExperiment`. Do or obtain the annotation data for the 
         assay in order to have the genomic coordinates ranges for all the features \n")
  
  # Keep the complete genomic coordinates in a GRanges object
  rse.ranges <- rowRanges(rangedsummarizedexperiment)
  
  # Find which genes in the RNA-seq assay are within our genomic coordinates
  genes.coords <- subsetByOverlaps(rse.ranges, gen.coords, type = "within")
  # Obtain their names
  genes.names <- names(genes.coords)
  genes.names
  
}

# Subset of Methylation CpGs within the specific genomic coordinates for a SummarizedExperiment object
# Must specify the name of the column where the genomic coordinate is

norangedCoords <- function(summarizedexperiment, gen.coords, col.name){
  
  # Check if there is rowData (annotation)
  if (length(rowData(summarizedexperiment)) == 0)
    stop("There is no rowData in the given `SummarizedExperiment`. Do or obtain the annotation data for the Methylation 
         assay in order to have the genomic coordinates for all the CpGs. \n")
  
  # Check Genomic Coordinates column from rowData is numeric or integer, if not convert
  if (! class(rowData(summarizedexperiment)[[col.name]]) %in% c("numeric","integer"))
    rowData(summarizedexperiment)[[col.name]] <- as.numeric(rowData(summarizedexperiment)[[col.name]])
  
  # Find which CpGs in the Methylation assay are within our genomic coordinates, and obtain their names
  cpgs.names <- rownames(rowData(summarizedexperiment)[rowData(summarizedexperiment)[[col.name]] >= start(gen.coords) & rowData(summarizedexperiment)[[col.name]] <= end(gen.coords), ])
  cpgs.names
  
}  


### NOTES ###

# Una función por cada omica/assay para hacer el subset 
# Qué hacer con la anotación? Si el assay no tiene rowData, hay que hacer la anotación, porque sino no se puede hacer
# el subset. Es mejor decirle al usuario que lo haga él o implementarlo en la función? Lo puedo poner en el ejemplo.
  



#################### OLD VERSION FUNCTIONS ##########################

# type.omics must be a character vector indicating the type of omics data from each assay from the multiassayexperiment
# in the CORRECT ORDER, e.g. type.omics = c("methylation","rnaseq")

genCoords_old <- function(multiassayexperiment, gen.coords, type.omics, col.name){
  
  # Check that the input is a MultiAssayExperiment
  if (!class(multiassayexperiment) == "MultiAssayExperiment")
    stop("Input must be a 'MultiAssayExperiment' object \n")
  
  # Check that gen.coords is a GRanges object
  if (!class(gen.coords) == "GRanges")
    stop("'gen.coords' must be a 'GRanges' object \n")
  
  ### Check that gen.coords and type.omics are provided.
  ### Check type.assay
  all.omics <- c("methylation", "proteomic", "genexp", "rnaseq")
  omics.match <- charmatch(type.omics, all.omics)
  
  subset.names = list()
  
  for (assay in 1:length(type.omics)) {
    
    ### Check which type of omics data is each assay. HOW TO DO THIS WITHOUT type.omics????
    # Methylation
    if (omics.match[assay] == 1) {
      # Obtain cpgs names
      subset.names[[names(experiments(multiassayexperiment)[assay])]] <- methCoords(multiassayexperiment[[assay]], gen.coords, col.name = "Genomic_Coordinate")
    }
    # RNA-seq
    else if (omics.match[assay] == 4) {
      # Obtain gene names
      subset.names[[names(experiments(multiassayexperiment)[assay])]] <- rnaseqCoords(multiassayexperiment[[assay]], gen.coords)      
    }
  }
  
  # MultiAssayExperiment subset using feature names obtained above with our genomic coordinate
  mae.subset <- subsetByRow(multiassayexperiment, subset.names)
  mae.subset
  
}

# Subset of RNA-seq genes within the specific genomic coordinates for RangedSummarizedExperiment object

rnaseqCoords_old <- function(rangedsummarizedexperiment, gen.coords){
  
  # Check if there is rowData (annotation)
  if (length(rowRanges(rangedsummarizedexperiment)) == 0)
    stop("There is no rowRanges data in the given `RangedSummarizedExperiment`. Do or obtain the annotation data for the RNA-seq 
         assay in order to have the genomic coordinates ranges for all the genes \n")
  
  # Keep the complete genomic coordinates in a GRanges object
  rse.ranges <- rowRanges(rangedsummarizedexperiment)
  
  # Find which genes in the RNA-seq assay are within our genomic coordinates
  genes.coords <- subsetByOverlaps(rse.ranges, gen.coords, type = "within")
  # Obtain their names
  genes.names <- names(genes.coords)
  genes.names
  
}

# Subset of Methylation CpGs within the specific genomic coordinates for a SummarizedExperiment object

methCoords_old <- function(summarizedexperiment, gen.coords, col.name){
  
  # Check if there is rowData (annotation)
  if (length(rowData(summarizedexperiment)) == 0)
    stop("There is no rowData in the given `SummarizedExperiment`. Do or obtain the annotation data for the Methylation 
         assay in order to have the genomic coordinates for all the CpGs. \n")
  
  # Check Genomic Coordinates column from rowData is numeric or integer, if not convert
  if (! class(rowData(summarizedexperiment)[[col.name]]) %in% c("numeric","integer"))
    rowData(summarizedexperiment)[[col.name]] <- as.numeric(rowData(summarizedexperiment)[[col.name]])
  
  # Find which CpGs in the Methylation assay are within our genomic coordinates, and obtain their names
  cpgs.names <- rownames(rowData(summarizedexperiment)[rowData(summarizedexperiment)[[col.name]] >= start(gen.coords) & rowData(summarizedexperiment)[[col.name]] <= end(gen.coords), ])
  cpgs.names
  
}  

  
#' Subset a MultiAssayExperiment object by specific genomic coordinates
#'
#' @description Keeps, in every assay of a \code{MultiAssayExperiment}, only the
#'   features that fall inside a genomic region of interest. Assays stored as a
#'   \code{RangedSummarizedExperiment} are filtered on their \code{rowRanges}
#'   (features contained \emph{within} the region); the assay indicated by
#'   \code{meth.index} -- typically methylation, held in a plain
#'   \code{SummarizedExperiment} without ranges -- is filtered on a numeric
#'   position column of its \code{rowData}. The resulting feature names are then
#'   used to subset the whole object row-wise.
#'
#' @param multiassayexperiment a \code{MultiAssayExperiment} object whose assays
#'   carry annotation (\code{rowRanges}, or \code{rowData} for the assay pointed
#'   at by \code{meth.index}).
#' @param gen.coords must be a GRanges object
#' @param meth.index is the index where the methylation data (or a SummarizedExperiment non Ranged) is located inside the mae.
#' @param col.name column name
#' @return A \code{MultiAssayExperiment} with the same assays as
#'   \code{multiassayexperiment}, each restricted to the features lying inside
#'   \code{gen.coords}.
#' @examples
#' if (requireNamespace("MultiAssayExperiment", quietly = TRUE) &&
#'     requireNamespace("SummarizedExperiment", quietly = TRUE) &&
#'     requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   library(MultiAssayExperiment)
#'   library(SummarizedExperiment)
#'   library(GenomicRanges)
#'
#'   data(cardiovascular)
#'   sel <- rownames(X2)[1:20]
#'
#'   ## a non-ranged assay (methylation): positions live in rowData
#'   m1 <- t(as.matrix(X1[rownames(X1) %in% sel, 1:8, drop = FALSE]))
#'   se1 <- SummarizedExperiment(
#'       assays  = list(beta = m1),
#'       rowData = DataFrame(Genomic_Coordinate = seq(1000, 8000,
#'                                                    length.out = nrow(m1)),
#'                           chr = "chr1",
#'                           row.names = rownames(m1)))
#'
#'   ## a ranged assay: positions live in rowRanges
#'   m2 <- t(as.matrix(X3[rownames(X3) %in% sel, , drop = FALSE]))
#'   gr <- GRanges(seqnames = "chr1",
#'                 ranges = IRanges(start = c(1200, 2200, 3200, 9000),
#'                                  width = 50))
#'   names(gr) <- rownames(m2)
#'   rse2 <- SummarizedExperiment(assays = list(expr = m2), rowRanges = gr)
#'
#'   mae <- MultiAssayExperiment(
#'       ExperimentList(list(methylation = se1, cells = rse2)))
#'
#'   region <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1000,
#'                                                         end = 5000))
#'   sub <- genCoords(mae, region, meth.index = 1,
#'                    col.name = "Genomic_Coordinate")
#'   dim(sub[[1]])   # methylation restricted to the region
#'   dim(sub[[2]])   # cells restricted to the region
#' }
#' @export

genCoords <- function(multiassayexperiment, gen.coords, meth.index, col.name){

  # Check that the input is a MultiAssayExperiment
  if (!inherits(multiassayexperiment, "MultiAssayExperiment"))
    stop("Input must be a 'MultiAssayExperiment' object \n")

  # Check that gen.coords is a GRanges object
  if (!inherits(gen.coords, "GRanges"))
    stop("'gen.coords' must be a 'GRanges' object \n")

  subset.names <- list()

  for (assay in seq_along(multiassayexperiment)) {

    # Check which type of omics data is each assay
    # No RangedSummarizedExperiment (e.g. Methylation)
    if (assay == meth.index) {
      # Obtain CpGs names
      subset.names[[names(MultiAssayExperiment::experiments(multiassayexperiment)[assay])]] <- norangedCoords(multiassayexperiment[[assay]], gen.coords, col.name = "Genomic_Coordinate")
    }
    # RangedSummarizedExperiment (e.g. RNA-seq, miRNA, proteomics, etc)
    else {
      # Obtain feature names
      subset.names[[names(MultiAssayExperiment::experiments(multiassayexperiment)[assay])]] <- rangedCoords(multiassayexperiment[[assay]], gen.coords)
    }
  }

  # MultiAssayExperiment subset using feature names obtained above with our genomic coordinate
  mae.subset <- MultiAssayExperiment::subsetByRow(multiassayexperiment, subset.names)
  mae.subset

}

# Subset of features within the specific genomic coordinates for a RangedSummarizedExperiment object

rangedCoords <- function(rangedsummarizedexperiment, gen.coords){

  # Check if there is rowData (annotation)
  if (length(SummarizedExperiment::rowRanges(rangedsummarizedexperiment)) == 0)
    stop("There is no rowRanges data in the given `RangedSummarizedExperiment`. Do or obtain the annotation data for the
         assay in order to have the genomic coordinates ranges for all the features \n")

  # Keep the complete genomic coordinates in a GRanges object
  rse.ranges <- SummarizedExperiment::rowRanges(rangedsummarizedexperiment)

  # Find which genes in the RNA-seq assay are within our genomic coordinates
  genes.coords <- IRanges::subsetByOverlaps(rse.ranges, gen.coords, type = "within")
  # Obtain their names
  genes.names <- names(genes.coords)
  genes.names

}

# Subset of Methylation CpGs within the specific genomic coordinates for a SummarizedExperiment object
# Must specify the name of the column where the genomic coordinate is

norangedCoords <- function(summarizedexperiment, gen.coords, col.name){

  # Check if there is rowData (annotation)
  if (length(SummarizedExperiment::rowData(summarizedexperiment)) == 0)
    stop("There is no rowData in the given `SummarizedExperiment`. Do or obtain the annotation data for the Methylation
         assay in order to have the genomic coordinates for all the CpGs. \n")

  # Check Genomic Coordinates column from rowData is numeric or integer, if not convert
  if (! class(SummarizedExperiment::rowData(summarizedexperiment)[[col.name]]) %in% c("numeric","integer"))
    SummarizedExperiment::rowData(summarizedexperiment)[[col.name]] <- as.numeric(SummarizedExperiment::rowData(summarizedexperiment)[[col.name]])

  # Find which CpGs in the Methylation assay are within our genomic coordinates, and obtain their names
  cpgs.names <- rownames(SummarizedExperiment::rowData(summarizedexperiment)[SummarizedExperiment::rowData(summarizedexperiment)[[col.name]] >= BiocGenerics::start(gen.coords) & SummarizedExperiment::rowData(summarizedexperiment)[[col.name]] <= BiocGenerics::end(gen.coords), ])
  cpgs.names

}


### NOTES ###

# Una función por cada omica/assay para hacer el subset
# Qué hacer con la anotación? Si el assay no tiene rowData, hay que hacer la anotación, porque sino no se puede hacer
# el subset. Es mejor decirle al usuario que lo haga él o implementarlo en la función? Lo puedo poner en el ejemplo.



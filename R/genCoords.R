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
#' @param gen.coords a \code{GRanges} object holding \emph{one} region. Several
#'   ranges are an error: pass them one at a time.
#' @param meth.index is the index where the methylation data (or a SummarizedExperiment non Ranged) is located inside the mae.
#' @param col.name name of the \code{rowData} column of the assay pointed at by
#'   \code{meth.index} that holds the feature position, in the same coordinate
#'   system as \code{gen.coords}.
#' @param seqname.col name of the \code{rowData} column of that same assay that
#'   holds the chromosome, so that features at the same position on another
#'   chromosome are dropped. The default \code{NULL} looks for a column called
#'   (case-insensitively) \code{seqnames}, \code{seqname}, \code{chr},
#'   \code{chrom} or \code{chromosome}; when none is found the position filter
#'   is applied on its own and a warning is issued.
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
#'   ## a non-ranged assay (methylation): positions live in rowData, and the
#'   ## chromosome is read from the 'chr' column without having to name it
#'   m1 <- t(as.matrix(X1[rownames(X1) %in% sel, 1:8, drop = FALSE]))
#'   se1 <- SummarizedExperiment(
#'       assays  = list(beta = m1),
#'       rowData = DataFrame(position = seq(1000, 8000,
#'                                          length.out = nrow(m1)),
#'                           chr = rep(c("chr1", "chr2"),
#'                                     length.out = nrow(m1)),
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
#'   sub <- genCoords(mae, region, meth.index = 1, col.name = "position")
#'   dim(sub[[1]])   # methylation, region and chromosome
#'   dim(sub[[2]])   # cells restricted to the region
#' }
#' @export

genCoords <- function(multiassayexperiment, gen.coords, meth.index, col.name,
                      seqname.col = NULL){

  # Check that the input is a MultiAssayExperiment
  if (!inherits(multiassayexperiment, "MultiAssayExperiment"))
    stop("Input must be a 'MultiAssayExperiment' object \n")

  # Check that gen.coords is a GRanges object
  if (!inherits(gen.coords, "GRanges"))
    stop("'gen.coords' must be a 'GRanges' object \n")

  # A single region: start() / end() would recycle silently over several ranges
  if (length(gen.coords) != 1L)
    stop("'gen.coords' must hold exactly one range, got ", length(gen.coords),
         ". Subset it, or call 'genCoords' once per range. \n", call. = FALSE)

  subset.names <- list()

  for (assay in seq_along(multiassayexperiment)) {

    # Check which type of omics data is each assay
    # No RangedSummarizedExperiment (e.g. Methylation)
    if (assay == meth.index) {
      if (missing(col.name))
        stop("'col.name' is required: it names the rowData column holding the ",
             "position of the features of the assay given by 'meth.index' \n",
             call. = FALSE)
      # Obtain CpGs names
      subset.names[[names(MultiAssayExperiment::experiments(multiassayexperiment)[assay])]] <- norangedCoords(multiassayexperiment[[assay]], gen.coords, col.name = col.name, seqname.col = seqname.col)
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

norangedCoords <- function(summarizedexperiment, gen.coords, col.name,
                           seqname.col = NULL){

  # Check if there is rowData (annotation)
  row.data <- SummarizedExperiment::rowData(summarizedexperiment)
  if (length(row.data) == 0)
    stop("There is no rowData in the given `SummarizedExperiment`. Do or obtain the annotation data for the Methylation
         assay in order to have the genomic coordinates for all the CpGs. \n")

  if (is.null(col.name) || !(col.name %in% colnames(row.data)))
    stop("'col.name' must name a column of the rowData of the assay given by ",
         "'meth.index'; '", col.name, "' is not one of: ",
         paste(colnames(row.data), collapse = ", "), " \n", call. = FALSE)

  # Check Genomic Coordinates column from rowData is numeric or integer, if not convert
  positions <- row.data[[col.name]]
  if (!is.numeric(positions))
    positions <- as.numeric(positions)

  # Find which CpGs in the Methylation assay are within our genomic coordinates
  keep <- positions >= BiocGenerics::start(gen.coords) &
          positions <= BiocGenerics::end(gen.coords)
  keep[is.na(keep)] <- FALSE

  # Position alone does not identify a locus: drop features sitting at the same
  # coordinate on a different chromosome.
  if (is.null(seqname.col))
    seqname.col <- .seqnameColumn(colnames(row.data))
  if (is.null(seqname.col)) {
    warning("no chromosome column found in the rowData of the non-ranged ",
            "assay, so features are selected on position alone. Give ",
            "'seqname.col' to filter on the chromosome as well.",
            call. = FALSE)
  } else {
    if (!(seqname.col %in% colnames(row.data)))
      stop("'seqname.col' must name a column of the rowData of the assay ",
           "given by 'meth.index'; '", seqname.col, "' is not one of: ",
           paste(colnames(row.data), collapse = ", "), " \n", call. = FALSE)
    keep <- keep & as.character(row.data[[seqname.col]]) ==
                   as.character(SummarizedExperiment::seqnames(gen.coords))
    keep[is.na(keep)] <- FALSE
  }

  # Obtain their names
  cpgs.names <- rownames(row.data)[keep]
  cpgs.names

}

# The chromosome column of a rowData(), looked up by the usual spellings.
# Returns NULL when none of them is present.

.seqnameColumn <- function(nms){
  candidates <- c("seqnames", "seqname", "chr", "chrom", "chromosome")
  hit <- match(candidates, tolower(nms))
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) NULL else nms[hit[1L]]
}


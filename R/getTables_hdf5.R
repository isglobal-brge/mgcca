#' Splits a MultiAssayExperiment in a list of tables
#'
#' @param multiassayexperiment MultiAssayExperiment
#' @param filename string with file name
#' @param overwriteFile boolean, if true , file is overwritten
#' @param overwriteDataset boolean, if true , dataset is overwritten
#' @return List of tables
#' @export
getTables_hdf5 <- function(multiassayexperiment, filename, overwriteFile = FALSE, overwriteDataset = FALSE){

  # Check if route exists
  if( !dir.exists(basename(dirname(filename))) ) {
    stop("Destination path does not exists")
  }
  # Check that the input is a MultiAssayExperiment
  if (!class(multiassayexperiment) == "MultiAssayExperiment")
    stop("Input must be a 'MultiAssayExperiment' object \n")

  ## Aquí s'haurà de posar una funció per a convertir de muytiassayExperiment
  ## cap a taules hdf5 --> Potser només cal utilitzar les funcions de
  ## crear fitxer i crear matrius ..... amb això potser n'hi hauria prou.....
  ## però potser caldria repassar les opcions de crear matrius amb dimnames i les opcions
  ## de restaurar els noms dels dimnames.
  ##

  #..# --  CODI ORIGINAL
  #..# tables.list = list()
  #..# for (assay in 1:length(multiassayexperiment)) {
  #..#   matrix.add <- as.matrix(assays(multiassayexperiment)[[assay]])
  #..#   tables.list[[names(assays(multiassayexperiment)[assay])]] <- t(matrix.add)
  #..# }
  #..# --  FI CODI ORIGINAL
  tables.list <-  list()
  for (assay in 1:length(experiments(multiassayexperiment))) {
    if(assay == 1){
      if( file.exists(filename) && overwriteFile == F) {
        bdAdd_hdf5_matrix(as.matrix(assays(multiassayexperiment)[[assay]]), filename, "MGCCA_IN", names(assays(multiassayexperiment)[assay]), force = overwriteDataset );
        tables.list[[assay]] <- names(assays(multiassayexperiment)[assay])
      } else {
        bdCreate_hdf5_matrix_file(filename, as.matrix(assays(multiassayexperiment)[[assay]]), "MGCCA_IN", names(assays(multiassayexperiment)[assay]), force = overwriteFile )
        tables.list[[assay]] <- names(assays(multiassayexperiment)[assay])
      }
      ##..## rownames(assays(multiassayexperiment)[[1]])
      ##..## colnames(assays(multiassayexperiment)[[1]])
    } else {
      bdAdd_hdf5_matrix(as.matrix(assays(multiassayexperiment)[[assay]]), filename, "MGCCA_IN", names(assays(multiassayexperiment)[assay]), force = overwriteDataset );
      tables.list[[assay]] <- names(assays(multiassayexperiment)[assay])}
  }

  #..# class(tables.list) <- "ListMAE"

  # Aquí potser estaria bé retornar els noms de les taules creades...
  # es a dir els noms dels datasets que s'han creat ??
  tables.list

}

print.listMAE <- function(listMAE) {
  print(paste("Object of class ", class(listMAE), sep = ""))
  print(paste(length(listMAE), " assays in the MultiAssayExperiment list:", sep = ""))
  print(names(listMAE))
  }



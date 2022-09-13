getScores_hdf5 <- function(filename, XX, As, threads) {

    bdapply_Function_hdf5(filename = filename,
                          group = "X", datasets = XX,
                          b_group = "scores/As", b_datasets = As,
                          outgroup = "scores/scores", func = "blockmult",transp_dataset = T,
                          force = T)
    scores <- bdgetDatasetsList_hdf5(filename = filename, group = "scores")
    return(scores)

}

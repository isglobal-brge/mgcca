getScores_hdf5 <- function(filename, XX, As, threads) {

    grp_score <- "FINAL_RESULTS/scores"
    bdapply_Function_hdf5(filename = filename,
                          group = "X", datasets = XX,
                          b_group = "scores/As", b_datasets = As,
                          outgroup = grp_score, func = "blockmult",transp_dataset = T,
                          force = T)
    scores <- bdgetDatasetsList_hdf5(filename = filename, group = grp_score)
    return(scores)

}

getScores_hdf5 <- function(i, dat, As){
    # scores <- dat[[i]]%*%As[[i]]


    #..# scores <- blockmult(dat[[i]],As[[i]], onmemory = T)
    #..# colnames(scores) <- paste0("comp", 1:ncol(scores))
    #..# scores

    bdapply_Function_hdf5(filename = filename,
                          group = "X",datasets = dat,
                          b_group = "As", b_datasets = As,
                          outgroup = "scores",func = "blockmult",
                          force = T)
    scores <- bdgetDatasetsList_hdf5(filename = filename, group = "scores")
    return(scores)


}

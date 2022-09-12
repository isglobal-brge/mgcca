getWeights_hdf5 <- function(A, XX, K, KX, threads){

    #..# KXA <- bdblockmult(bdblockmult(K[[i]],XX[[i]], onmemory = T),A[[i]], onmemory = T)
    bdapply_Function_hdf5(filename = filename,
                          group = "KX",datasets = KX,
                          b_group = "A", b_datasets = A,
                          outgroup = "KXA",func = "blockmult",
                          force = T)
    KXA <- bdgetDatasetsList_hdf5(filename = filename, group = "KXA")

    #..# vv <- 1/apply(KXA, 2, sd)
    # vv <- 1/apply(KXA, 2, sd)
    bdapply_Function_hdf5(filename = filename,
                          group = "KXA",datasets = KXA,
                          outgroup = "vv",func = "sdmean",
                          force = T)

    browser()
    #..#    vv <-  bdblockmult(cbind(rep(1,nrow(A[[i]]))),rbind(vv) , onmemory = T)
    #..#    As <- A[[i]]*vv
    #..#    rownames(As) <- colnames(XX[[i]])
    #..#    As


}

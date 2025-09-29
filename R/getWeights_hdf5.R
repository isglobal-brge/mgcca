getWeights_hdf5 <- function(filename, A, XX, K, KX, initialGroup, threads){

    bdapply_Function_hdf5(filename = filename,
                          group = "KX",datasets = KX,
                          b_group = "scores/A", b_datasets = A,
                          outgroup = "KXA", func = "blockmult",
                          force = T)
    KXA <- bdgetDatasetsList_hdf5(filename = filename, group = "KXA")

    bdapply_Function_hdf5(filename = filename,
                          group = "KXA",datasets = KXA,
                          outgroup = "vv",func = "sdmean",
                          force = T)

    vvl <- bdgetDatasetsList_hdf5(filename = filename, group = "vv", prefix = "sd")

    sapply(1:length(A), function(i) {
        Ad <- h5read(filename, paste0("scores/A/",A[i]))
        vv <- 1/(rhdf5::h5read(filename, paste0("vv/",vvl[i])))
        vv <-  bdblockmult(cbind(rep(1,nrow(Ad))), t(vv) , onmemory = T)
        As <-  Ad*vv
        rownames(As) <- rhdf5::h5read(filename, paste0(initialGroup,"/.",XX[i],"_dimnames/1"))[,1]
        bdAdd_hdf5_matrix(object = As, filename = filename,
                          group = "scores/As",dataset = XX[i],
                          force = T )
    })

    bdgetDatasetsList_hdf5(filename = filename, group = "scores/As")

}

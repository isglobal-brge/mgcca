getCor_hdf5 <- function(filename, Xgroup, x, Ygroup, y, byblocks, threads) {

    common <- BigDataStatMeth::bdgetDiagonal_hdf5(filename, group = "K", dataset =  x)
    common_elems <- which(common!=0)

    # Subset dataset if there are missing samples inside
    if( length(common_elems) != length(common) ) {
        xgroup_inter <- "tmp"
        ygroup_inter <- "tmp"

        bdsubset_hdf5_dataset(filename = filename, dataset_path = paste0( Ygroup, "/", y),
                              indices = common_elems,  select_rows = TRUE,
                              new_group = "tmp",  new_name = y, overwrite = TRUE)

        bdsubset_hdf5_dataset(filename = filename, dataset_path = paste0( Xgroup, "/", x),
                              indices = common_elems, select_rows = FALSE,
                              new_group = "tmp",  new_name = x, overwrite = TRUE)
    } else {
        xgroup_inter <- Xgroup
        ygroup_inter <- Ygroup
    }

    res <- bdCorr_hdf5( filename_x = filename, group_x = xgroup_inter, dataset_x = x, trans_x = TRUE,
            filename_y = filename, group_y = ygroup_inter, dataset_y = y, compute_pvalues = TRUE)

    bdmove_hdf5_dataset(filename,source_path = paste0(res$group, "/", res$correlation ),
                        dest_path =  paste0("FINAL_RESULTS/corsY/", x), overwrite = TRUE )
    bdWrite_hdf5_dimnames(filename = filename,
                       group = "FINAL_RESULTS/corsY/",
                       dataset = x,
                       rownames = t(getDimNames_hdf5(filename, Xgroup, x)$rownames),
                           colnames = paste0("comp", seq_len(res$n_variables_y)))

    bdmove_hdf5_dataset(filename,source_path = paste0(res$group, "/", res$pvalues ),
                        dest_path =  paste0("FINAL_RESULTS/pval.cor/", x), overwrite = TRUE )

    bdWrite_hdf5_dimnames(filename = filename,
                       group = "FINAL_RESULTS/pval.cor/",
                       dataset = x,
                       rownames = t(getDimNames_hdf5(filename, Xgroup, x)$rownames),
                       colnames = paste0("comp", seq_len(res$n_variables_y)))

    # Remove intermediate datasets
    if(xgroup_inter == "tmp") {
        BigDataStatMeth::bdRemove_hdf5_element(filename, paste0(xgroup_inter, "/", x ) )
        BigDataStatMeth::bdRemove_hdf5_element(filename, paste0(ygroup_inter, "/", y ) )
    }

}


#ifndef MGCCA_HPP
#define MGCCA_HPP

#include "mgcca.h"



/**
 * @brief Split HDF5 dataset paths into group name and dataset names.
 *
 * @param paths Vector of full dataset paths (e.g., "GROUP/DATASET").
 *
 * @return std::pair<std::string, std::vector<std::string>>
 *         - first  : group name (substring before the last '/')
 *         - second : dataset names (substrings after the last '/')
 *
 * @note If multiple paths have different groups, only the group of the
 *       last processed path is returned.
 *
 * @code
 * std::vector<std::string> v = {
 *   "MGCCA_IN/ACC_Methylation",
 *   "MGCCA_IN/ACC_RNASeq2GeneNorm-20160128"
 * };
 * auto res = splitPaths(v);
 * // res.first  == "MGCCA_IN"
 * // res.second == {"ACC_Methylation", "ACC_RNASeq2GeneNorm-20160128"}
 * @endcode
 */
 std::pair<std::string, std::vector<std::string>>
     splitPaths(const std::vector<std::string> &paths) {
         std::string group;
         std::vector<std::string> datasets;
         for (const auto &p : paths) {
             size_t pos = p.find_last_of('/');
             if (pos != std::string::npos) {
                 group = p.substr(0, pos);
                 datasets.push_back(p.substr(pos + 1));
             }
         }
         return {group, datasets};
     }



 /**
  * @brief MGCCA over HDF5 using BigDataStatMeth (C++ backend)
  * @details Compute Generalized CCA components on multi-omic tables stored
  * in an HDF5 file, operating block-wise via BigDataStatMeth. Results are
  * written under the group "MGCCA_OUT" in the same file.
  *
  * @param filename Path to the HDF5 file.
  * @param group Input HDF5 group that contains the tables.
  * @param datasets Dataset names inside @p group.
  * @param nfac Number of shared components to extract.
  * @param scale 0 = none; 1 = center; 2 = center+standardize.
  * @param pval Power parameter p for the GCCA objective.
  * @param scores 0/1; if 1, write per-table scores to HDF5.
  * @param method GCCA scheme (e.g., "horst", "centroid", "factorial").
  * @param lambda Ridge/shrinkage penalty (applied per table if scalar).
  * @param mc_cores Number of CPU cores/threads to use.
  *
  * @return void. Primary outputs are written to HDF5.
  *
  * @note C++ identifiers cannot contain dots. The R-side argument
  * "mccores" maps to @p mc_cores here.
  *
  * @par References
  * van de Velden, M., & Takane, Y. (2012). Generalized CCA with missing
  * values for multi-block data integration.
*/

//' @title MGCCA over HDF5 using BigDataStatMeth (C++ backend)
//' @description Compute Generalized CCA components on multi-omic tables
//'   stored in an HDF5 file, operating block-wise via BigDataStatMeth.
//'   Results are written under the group \code{MGCCA_OUT} in the same file.
//'
//' @param filename String. Path to the HDF5 file.
//' @param group String. Input HDF5 group that contains the tables.
//' @param datasets Character vector. Dataset names inside \code{group}.
//' @param nfac Integer. Number of shared components to extract.
//' @param scale Integer. 0 = none, 1 = center, 2 = center+standardize.
//' @param pval Integer. Power parameter \code{p} for the GCCA objective.
//' @param scores Integer (0/1). If 1, write per-table scores to HDF5.
//' @param method String. GCCA scheme (e.g., "horst", "centroid", "factorial").
//' @param lambda Numeric. Ridge/shrinkage penalty (applied per table if scalar).
//' @param mccores Integer. Number of CPU cores/threads to use.
//'
//' @details Inputs are read block-wise from HDF5; large intermediates are
//'   avoided. Outputs include \code{Y} (shared components), \code{corsY},
//'   \code{scores/<table>} and AVE metrics inside \code{MGCCA_OUT}.
//'
//' @return Invisibly returns \code{TRUE}; primary outputs are written to HDF5.
//'
//' @references van de Velden, M. & Takane, Y. (2012). Generalized CCA with
//'   missing values for multi-block data integration.
//'
//' @note C++ identifiers cannot contain dots. Implement the argument as
//'   \code{mc_cores} in C++ and expose/rename to \code{mccores} on the R side
//'   if you need the dotted name in R.
//'
//' @export
// [[Rcpp::export]]
void mgcca_rcpp( std::string filename, std::string group,
                 std::vector<std::string> datasets, int nfac, int scale,
                 int pval, int scores, std::string method, Rcpp::Nullable<std::vector<double>> lambda,
                 Rcpp::Nullable<int> mccores)
{


    hdf5File* dsFile = nullptr;
    hdf5Dataset* dstmp = nullptr;

    try {

        std::vector<double> dlambda;

        // Define valid inversion types
        Rcpp::IntegerVector inv_types = {1, 2};
        inv_types.names() = Rcpp::CharacterVector({"solve", "penalized"});

        int idx = inv_types.findName(method);
        if (idx == -1) stop("method should be 'solve' or 'penalized'");

        int n = datasets.size();


        // Fast size without copying:
        int nlambda = lambda.isNull() ? 0 : Rf_length(lambda.get());

        if (inv_types[idx] == 2 && (lambda.isNull() || nlambda != n )) {
            stop(lambda.isNull() ? "penalized method requires lambda parameter" :
                     "lambda must be a vector of length equal to the number of tables");
        }

        // input vector with full HDF5 paths
        std::vector<std::string> paths = datasets; // or rename your input to 'paths'

        if (group.empty()) {
            auto [grp, ds] = splitPaths(paths);
            group = std::move(grp);
            datasets = std::move(ds);
        }
        //
        // if( group == "" ) {
        //     auto res = splitPaths(datasets);
        //     group = res.first;      // group
        //     datasets = res.second;  // datasets
        //
        //     Rcpp::Rcout<<"Groups: "<< group<<"\n";
        //     for (auto &d : datasets) std::cout << d << ' ';
        //     std::cout << std::endl;
        // }


        // Create FINAL RESULTS folder
        hdf5Group* dsGroup = new hdf5Group(filename, "FINAL_RESULTS");
        delete dsGroup; dsGroup = nullptr;

        // if (scale) {
        // #..# x <- lapply(x, scale)
        //
        //     // RcppApplyFunctionHdf5(filename, )
        //     //
        //     //
        //     // ("MGCCA_IN/ACC_Methylation", "MGCCA_IN/ACC_RNASeq2GeneNorm-20160128")
        //     //     datasets = c("ACC_Methylation-20160128", "ACC_RNASeq2GeneNorm-20160128")
        //     //
        //     // RcppApplyFunctionHdf5( std::string filename,
        //     //                        std::string group,
        //     //                        Rcpp::StringVector datasets,
        //     //                        std::string outgroup,
        //     //                        std::string func,
        //     //                        Rcpp::Nullable<Rcpp::CharacterVector> b_group = R_NilValue,
        //     //                        Rcpp::Nullable<Rcpp::StringVector> b_datasets = R_NilValue,
        //     //                        Rcpp::Nullable<bool> overwrite = false,
        //     //                        Rcpp::Nullable<bool> transp_dataset = false,
        //     //                        Rcpp::Nullable<bool> transp_bdataset = false,
        //     //                        Rcpp::Nullable<bool> fullMatrix = false,
        //     //                        Rcpp::Nullable<bool> byrows = false,
        //     //                        Rcpp::Nullable<int> threads = R_NilValue)
        //     //
        //
        //     // x <- lapply(datasets, bdNormalize_hdf5, filename = filename, group = group, bcenter=TRUE, bscale=TRUE, force = TRUE, byrows = TRUE )
        // }

    } catch(std::exception& ex) {
        Rcout<< "c++ exception getQRbyBlocks_rcpp: "<< ex.what() << "\n";
        return void();
    }

    return void();

}





//  ----------------------------------------------------------------------------
//                      FUNCIÓ MGCCA A IMPLEMENTAR !!!
//  ----------------------------------------------------------------------------


//
//  mgcca_hdf5 <- function(x, filename, group, datasets, nfac=2, scale=TRUE, pval=TRUE, scores=FALSE,
//                            method="solve", lambda, mc.cores=1, ...)
//  {
//
// # #' @importFrom rfunctions geninv cgls
//         inv.type <- c("solve", "penalized")
//             inv.method <- charmatch(method, inv.type, nomatch = 0)
//             if (inv.method == 0)
//                 stop("method should be 'solve' or 'penalized' \n")
//
//                 n <- length(datasets) # number of tables
//
//             if (inv.method == 2) {
//                 if (missing(lambda)) {
//                     stop("penalized method requires lambda parameter \n")
//                 } else {
//                     if(length(lambda)!=n) {
//                         stop("lambda must be a vector of length equal to the number of tables \n")
//                     }
//                 }
//             }
//
// # Create FINAL RESULTS folder
//             BigDataStatMeth::bdCreateGroup_hdf5(filename, "FINAL_RESULTS" )
//
//                 if (scale) {
// #..# x <- lapply(x, scale)
//                     x <- lapply(datasets, bdNormalize_hdf5, filename = filename, group = group, bcenter=TRUE, bscale=TRUE, force = TRUE, byrows = TRUE )
//                 }
//
//                 currentdatasets <- sapply(datasets, function(d, g) {
//                     if(scale) {
//                         daasetname <- paste("NORMALIZED", g, d, sep = "/")
//                     } else {
//                         daasetname <- paste( g, d, sep = "/")
//                     }
//                 }, g = group)
//
//                     distancia <- regexpr("\\/[^\\/]*$", currentdatasets[1])[[1]]
//                 ngroup <-  substr(currentdatasets[1],1,distancia-1)
//
//                     if(length(currentdatasets)<=1)
//                         stop("we need more than one dataset to perform mgcca analysis")
//
//                         ns <- sapply( currentdatasets, function(dataset, file) { return( bdgetDim_hdf5(file, dataset)[2] ) }, file = filename )
//
//
//                         if(max(ns)==min(ns)) { # check whether there are missing individuals
//                             rn <- Reduce('union', sapply(datasets, getRowNames_hdf5, filename = filename, group = group))
//                         } else {
//                             rn <- sort(Reduce('union', sapply(datasets, getRowNames_hdf5, filename = filename, group = group)))
//                         }
//
//                         m <- length(rn)  # get the maximum number of individuals
//
//                         mclapply( datasets, getK_hdf5, ids=rn, m=m, mc.cores=mc.cores, filename = filename, group = group, ngroup = ngroup )
//                             X <- bdgetDatasetsList_hdf5(filename = filename, group = "X")
//                             K <- bdgetDatasetsList_hdf5(filename = filename, group = "K")
//                             p <-  sapply( paste0( "X/",X), function(el, filename){
//                                 res <- BigDataStatMeth::bdgetDim_hdf5(filename, el)
//                                 return(res[2]) # number of variables per table
//                             }, filename = filename )
//
//                             numvars <- min(p) # minimum number of variables
//
//                         getXKX_hdf5(filename, X, K, inv.method, lambda=lambda, scores, mc.cores=mc.cores)
//
//                             XKX <- bdgetDatasetsList_hdf5(filename = filename, group = "XKX")
//                             Mi <- solution_hdf5( filename = filename, X = X, XKX = XKX, mc.cores)
//
//                             bdReduce_matrix_hdf5(filename = filename, group = Mi, reducefunction = "+", outgroup = "FinalRes", outdataset = "M")
//
//                             bdReduce_matrix_hdf5(filename = filename, group = "K", reducefunction = "+", outgroup = "FinalRes", outdataset = "Ksum")
//
//                             Ksum <- bdgetDiagonal_hdf5(filename, "FinalRes", "Ksum")
//                             bdAdd_hdf5_matrix(as.matrix(Ksum^-0.5), filename, "FinalRes", "Ksum05")
//                             mult_wXw_hdf5(filename, "FinalRes","M", "FinalRes", "Ksum05")
//
//                             bdSVD_hdf5( filename,group = "FinalRes", dataset = "MKsum05",bcenter = F,bscale = F)
//
//                             Yast <- Re( rhdf5::h5read(filename,"SVD/MKsum05/u",))[,1:nfac]
//
//                         Y <- sqrt(n)* bdwproduct(Yast, Ksum^-0.5, "wX")
//                             colnames(Y) <- paste0("comp", 1:ncol(Y))
//                             rownames(Y) <- rn
//
//                         sapply(datasets, function(dataset) {
//                             bdAdd_hdf5_matrix(object = Y, filename = filename, group = "FinalRes/Y",
//                                               dataset = paste0(dataset,".Y"), force = T)
//                         } )
//
//                             if (scores) {
//                                 Yd <- bdgetDatasetsList_hdf5(filename = filename, group = "FinalRes/Y")
//                                 KX <- bdgetDatasetsList_hdf5(filename = filename, group = "KX")
//                                 XK <- bdgetDatasetsList_hdf5(filename = filename, group = "XK")
//                                 A <- productXKY_hdf5(filename=filename, Y=Yd, XK=XK, XKX=XKX, mc.cores)
//                                 As <- getWeights_hdf5(filename, A=A, XX=X, K=K, KX=KX, initialGroup = group, mc.cores)
//                                 scores <- getScores_hdf5(filename, XX = X, As = As, mc.cores)
//                                 sapply(scores, function(sdataset) {
//                                     bdWriteDimnames_hdf5( filename, group = "FINAL_RESULTS/scores", dataset = sdataset,rownames = rn,colnames = colnames(Y))
//                                 } )
//                             } else {
//                                 scores <- NULL
//                             }
//
//                             corsY <- lapply( X, function(dataX) {
//                                 res <- getCor_hdf5(filename = filename, x = dataX, Xgroup = "X", Ygroup = "FinalRes/Y",
//                                                    y = Y, byblocks = T, threads = mc.cores)
//                                 return(res)
//                             })
//
//                                 sapply(1:length(corsY), function(i) {
//                                     bdAdd_hdf5_matrix( corsY[[i]] , filename, group = "FINAL_RESULTS/corsY",
//                                                        dataset = X[i], force = T)
//                                 })
//
// # Create FINAL RESULTS
//
//                                 if (pval) {
//                                     pval.cor <- mclapply(corsY, cor.test.p, n=m, mc.cores=mc.cores)
//                                     sapply(1:length(pval.cor), function(i) {
//                                         bdAdd_hdf5_matrix( pval.cor[[i]] , filename, group = "FINAL_RESULTS/pval",
//                                                            dataset = X[i], force = T)
//                                     })
//                                 } else {
//                                     pval.cor <- NULL
//                                 }
//
//                                 AVE_X <- lapply(corsY, function(x) apply(x^2, 2, mean))
//                                     outer <- matrix(unlist(AVE_X), nrow = nfac)
//                                     AVE_outer <- sapply(1:nfac, function(j, p) sum(p * outer[j,])/sum(p),
//                                                         p=p)
// # AVE_inner <- Re(eig$v[1:nfac])
//                                     AVE_inner <- Re( rhdf5::h5read(filename,"SVD/MKsum05/v",))[,1:nfac]
//
//                                 sapply(1:length(AVE_X), function(i) {
//                                     bdAdd_hdf5_matrix( as.matrix(AVE_X[[i]]) , filename, group = "FINAL_RESULTS/AVE/AVE_X",
//                                                        dataset = X[i], force = T)
//                                 })
//                                     sapply(1:length(AVE_outer), function(i) {
//                                         bdAdd_hdf5_matrix( AVE_outer[[i]] , filename, group = "FINAL_RESULTS/AVE/AVE_outer",
//                                                            dataset = X[i], force = T)
//                                     })
//
//                                     bdAdd_hdf5_matrix( AVE_inner , filename, group = "FINAL_RESULTS/AVE",
//                                                        dataset = "AVE_inner", force = T)
//
//                                     bdAdd_hdf5_matrix( Y , filename, group = "FINAL_RESULTS",
//                                                        dataset = "Y", force = T)
//  }



#endif // MGCCA_HPP

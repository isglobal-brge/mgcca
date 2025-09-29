#include "BigDataStatMeth.hpp"

/**
 * @file get_AVE.cpp
 * @brief Implementation of Average Variance Extracted (AVE) calculation for GCCA
 */

/**
 * @brief Calculate Average Variance Extracted (AVE) for multiple datasets in HDF5 file
 *
 * @details This function computes the Average Variance Extracted (AVE) for multiple datasets
 * stored in HDF5 format as part of Generalized Canonical Correlation Analysis (GCCA).
 * The function processes each dataset by reading data in blocks, computing squared values
 * and row-wise means, then combines results across datasets to compute both inner and
 * outer model AVE statistics.
 *
 * The function performs the following operations:
 * 1. Iterates through each specified dataset in the HDF5 group
 * 2. Reads data in memory-efficient blocks to handle large datasets
 * 3. Computes squared values and row-wise means for each block
 * 4. Stores individual dataset AVE results in /FINAL_RESULTS/AVE/AVE_X/[dataset_name]
 * 5. Computes weighted outer model AVE using dataset sizes as weights
 * 6. Stores outer model AVE in /FINAL_RESULTS/AVE/AVE_outer_model
 * 7. Moves inner model AVE from /EIGEN/MKsum05/values to /FINAL_RESULTS/AVE/AVE_inner_model
 *
 * @param filename Path to the HDF5 file containing the datasets
 * @param group Group path within HDF5 file where datasets are located (e.g., "/FINAL_RESULTS/corsY")
 * @param datasets Vector of dataset names to process within the specified group
 * @param nfac Number of factors/components in the GCCA model
 * @param p Vector containing the number of variables for each dataset (used as weights)
 * @param wsize Optional block size for memory-efficient processing. If not specified,
 *              default block size of 1000 is used
 *
 * @return Rcpp::List Currently returns an empty list (implementation incomplete)
 *
 * @throws H5::FileIException If there are issues accessing the HDF5 file
 * @throws H5::DataSetIException If there are issues with HDF5 dataset operations
 * @throws std::exception For any other standard C++ exceptions
 *
 * @note The function creates new datasets in the HDF5 file under /FINAL_RESULTS/AVE/
 * @note Memory usage is controlled through block-wise processing of large datasets
 * @note The function assumes datasets contain numeric data suitable for AVE calculation
 *
 * @see van de Velden, M., & Takane, Y. (2012). Generalized canonical correlation analysis
 *      with missing values. Computational Statistics, 27(4), 551-571.
 *
 * @example
 * ```cpp
 * std::vector<std::string> datasets = {"dataset1", "dataset2", "dataset3"};
 * std::vector<int> p = {100, 150, 200};
 * auto result = get_AVE("data.h5", "/FINAL_RESULTS/corsY", datasets, 3, p, 500);
 * ```
 */

//' Calculate Average Variance Extracted (AVE) for GCCA Multi-omic Data
//'
//' @description
//' Computes the Average Variance Extracted (AVE) statistics for multiple datasets
//' stored in HDF5 format as part of Generalized Canonical Correlation Analysis (GCCA)
//' with missing values. This function is designed for multi-omic data integration
//' and provides both inner and outer model AVE statistics.
//'
//' @details
//' The Average Variance Extracted (AVE) is a measure of the amount of variance
//' captured by the latent components relative to measurement error. This function
//' calculates AVE for multiple datasets (e.g., different omics layers) and provides:
//'
//' \itemize{
//'   \item Individual dataset AVE statistics stored in HDF5 format
//'   \item Weighted outer model AVE combining all datasets
//'   \item Inner model AVE from eigenvalue decomposition
//' }
//'
//' The function processes large datasets efficiently using block-wise reading
//' to minimize memory usage while maintaining computational accuracy.
//'
//' @param filename Character string specifying the path to the HDF5 file
//' @param group Character string specifying the group path within the HDF5 file
//'        where the datasets are located (e.g., "/FINAL_RESULTS/corsY")
//' @param datasets Character vector containing the names of datasets to process
//' @param nfac Integer specifying the number of factors/components in the GCCA model
//' @param p Integer vector containing the number of variables for each dataset.
//'        Used as weights in the outer model AVE calculation
//' @param wsize Integer specifying the block size for memory-efficient processing.
//'        Default is NULL, which uses an internal default of 1000
//'
//' @return
//' Currently returns an empty list. The actual results are stored directly in the
//' HDF5 file under the following locations:
//' \itemize{
//'   \item \code{/FINAL_RESULTS/AVE/AVE_X/[dataset_name]} - Individual dataset AVE
//'   \item \code{/FINAL_RESULTS/AVE/AVE_outer_model} - Weighted outer model AVE
//'   \item \code{/FINAL_RESULTS/AVE/AVE_inner_model} - Inner model AVE from eigenvalues
//' }
//'
//' @section HDF5 Structure:
//' The function expects the input datasets to be stored in the specified group
//' within the HDF5 file. Each dataset should contain numeric data organized
//' as a matrix where rows represent observations and columns represent variables.
//'
//' @section Memory Management:
//' Large datasets are processed in blocks to control memory usage. The block size
//' can be controlled via the \code{wsize} parameter. Smaller block sizes use less
//' memory but may be slower for very large datasets.
//'
//' @section Error Handling:
//' The function includes comprehensive error handling for HDF5 file operations
//' and will provide informative error messages if issues occur during processing.
//'
//' @references
//' van de Velden, M., & Takane, Y. (2012). Generalized canonical correlation analysis
//' with missing values. Computational Statistics, 27(4), 551-571.
//'
//' @examples
//' \dontrun{
//' # Example usage for multi-omic GCCA analysis
//' library(mgcca)
//'
//' # Define datasets and parameters
//' datasets <- c("genomics", "transcriptomics", "proteomics")
//' p_values <- c(1000, 2000, 500)  # Number of variables per dataset
//' n_factors <- 3
//'
//' # Calculate AVE statistics
//' result <- get_AVE(
//'   filename = "multiomics_data.h5",
//'   group = "/FINAL_RESULTS/corsY",
//'   datasets = datasets,
//'   nfac = n_factors,
//'   p = p_values,
//'   wsize = 1000
//' )
//'
//' # Results are stored in HDF5 file and can be accessed using:
//' # - /FINAL_RESULTS/AVE/AVE_X/ for individual dataset AVE
//' # - /FINAL_RESULTS/AVE/AVE_outer_model for combined outer model AVE
//' # - /FINAL_RESULTS/AVE/AVE_inner_model for inner model AVE
//' }
//'
//' @seealso
//' \code{\link{mgcca}}, \code{\link{gcca_missing}}
//'
//' @author Your Name
//' @export
 // [[Rcpp::export]]
Rcpp::List get_AVE(const std::string& filename,
                   const std::string& group,
                   const std::vector<std::string>& datasets,
                   const int nfac, std::vector<int>& p,
                   Rcpp::Nullable<int> wsize = R_NilValue)
{

    BigDataStatMeth::hdf5Dataset* dsA = nullptr;
    BigDataStatMeth::hdf5Dataset* dsOut = nullptr;

    Rcpp::List results(datasets.size());
    Rcpp::CharacterVector names(datasets.size());

    Eigen::MatrixXd outer = Eigen::MatrixXd::Zero(nfac,datasets.size());

    try {

        int block_size = 1000;
        hsize_t nrows, ncols;

        std::vector<hsize_t> stride = {1, 1},
                             block = {1, 1};


        // Iterate through each dataset
        for(size_t i = 0; i < datasets.size(); ++i) {

            std::vector<hsize_t>  offset = {0, 0},
                                  count = {0, 0};

            // Create hdf5Dataset object for current dataset
            dsA = new BigDataStatMeth::hdf5Dataset(filename, group, datasets[i], false);
            dsA->openDataset();

            hsize_t* dims_out = dsA->dim();

            // Get dataset dimensions
            nrows = dsA->nrows();
            ncols = dsA->ncols();

            // Initialize to store results
            Eigen::MatrixXd datamean = Eigen::MatrixXd::Zero(1,nrows);

            // block_size = BigDataStatMeth::get_block_size(wsize, dims_out[1], dims_out[0]);

            block_size = std::min( (int)nrows, block_size);
            count[1] = dims_out[1];
            count[0] = block_size;

            // Read data in blocks of 500 columns
            for(hsize_t i=0; (i <= floor(dims_out[0]/block_size)) || i==0; i++)
            {

                if( offset[0] + block_size <= dims_out[0] ) {
                    count[0] = block_size;
                }else {
                    count[0] = dims_out[0] - offset[0];
                }

                std::vector<double> vdA( count[0] * count[1] );
                dsA->readDatasetBlock( {offset[0], offset[1]}, {count[0], count[1]}, stride, block, vdA.data() );
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X (vdA.data(), count[0], count[1] );

                Eigen::VectorXd sqmean = X.array().square().rowwise().mean();

                datamean.block( 0, offset[0], 1, sqmean.size()) = sqmean.transpose();

                offset[0] = offset[0] + block_size;
            }

            // Set outer
            outer.col(i) = datamean.row(0).transpose();

            // Create output dataset
            dsOut = new BigDataStatMeth::hdf5Dataset(filename, "/FINAL_RESULTS/AVE/AVE_X", datasets[i], true);
            // dsOut->createDataset(nrows, 1, "real");
            dsOut->createDataset(nrows, 1, "real");

            dsOut->writeDataset(Rcpp::wrap(datamean));

            // Set colnames
            Rcpp::CharacterVector colnames(nfac);
            for(int i = 0; i < nfac; ++i) {
                colnames[i] = "comp" + std::to_string(i + 1);
            }

            BigDataStatMeth::hdf5Dims* dsdims = nullptr;
            dsdims = new BigDataStatMeth::hdf5Dims(dsOut);

            if( dsdims != nullptr) {
                Rcpp::StringVector svrownames(1);
                dsdims->writeDimnames( Rcpp::wrap(svrownames), colnames);
                delete dsdims; dsdims = nullptr;
            }

            // Clean up current datasets
            delete dsA; dsA = nullptr;
            delete dsOut; dsOut = nullptr;

        }

        Eigen::Map<Eigen::VectorXi> p_map(p.data(), p.size());
        Eigen::VectorXd p_double = p_map.cast<double>();

        double sum_p = p_double.sum();
        Eigen::VectorXd AVE_outer = (outer.array().colwise() * p_double.array()).rowwise().sum() / sum_p;

        dsOut = new BigDataStatMeth::hdf5Dataset(filename, "/FINAL_RESULTS/AVE", "AVE_outer_model", true);
        dsOut->createDataset(AVE_outer.size(), 1, "real");
        dsOut->writeDataset(Rcpp::wrap(AVE_outer));

        delete dsOut; dsOut = nullptr;

        dsA = new BigDataStatMeth::hdf5Dataset(filename, "/EIGEN/MKsum05", "values", false);
        dsA->openDataset();
        dsA->moveDataset("/FINAL_RESULTS/AVE/AVE_inner_model", true);

        // Clean up current datasets
        delete dsA; dsA = nullptr;

    } catch(H5::FileIException& error) {
        checkClose_file(dsA, dsOut);
        Rf_error("c++ exception get_AVE (File IException)");
    } catch(H5::DataSetIException& error) {
        checkClose_file(dsA, dsOut);
        Rf_error("c++ exception get_AVE (DataSet IException)");
    } catch(std::exception& error) {
        checkClose_file(dsA, dsOut);
        Rf_error("c++ exception get_AVE function: %s", error.what());
    }

    return Rcpp::List(); // Should never reach here
}

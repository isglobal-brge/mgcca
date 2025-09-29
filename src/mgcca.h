#ifndef MGCCA_H
#define MGCCA_H

#include <Rcpp.h>
#include "BigDataStatMeth.hpp"


using namespace Rcpp;
using namespace BigDataStatMeth;

    std::pair<std::string, std::vector<std::string>> splitPaths(const std::vector<std::string> &paths);
    void mgcca_rcpp( std::string filename, std::string group, std::vector<std::string> datasets, int nfac, int scale,
                     int pval, int scores, std::string method, double lambda, int mccores);


#endif // MGCCA_H

##########################
# User-defined functions #
##########################

# Cpp functions
library(Rcpp)

suppressMessages(sourceCpp("DetectAndSolve_debugged.cpp"))

suppressMessages(
    sourceCpp(
        code = "// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
                #include <RcppArmadillo.h>
                #include <RcppEigen.h>

                // [[Rcpp::export]]
                SEXP armaMatMult(arma::mat A, arma::mat B){
                    arma::mat C = A * B;

                    return Rcpp::wrap(C);
                }

                // [[Rcpp::export]]
                SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
                    Eigen::MatrixXd C = A * B;

                    return Rcpp::wrap(C);
                }

                // [[Rcpp::export]]
                SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
                    Eigen::MatrixXd C = A * B;

                    return Rcpp::wrap(C);
        }",
        showOutput = FALSE,
        verbose = FALSE,
        echo = FALSE
    )
)

# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

# Standardize
Standardize <- function(M) {
    # Centralize
    M <- M - matrix(rep(colMeans(M), times = nrow(M)), nrow = nrow(M) , ncol = ncol(M), byrow = T)

    # Standardize
    M <- sweep(M, 2, sqrt(apply(M, 2, crossprod) / nrow(M)), "/")

    return(M)
}

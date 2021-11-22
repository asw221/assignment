
#include <vector>

#include <Rcpp.h>
#include <RcppEigen.h>

#include "assignment_problem.hpp"



/* "Assignment Problem" cost minimization
 *
 * Find minimum assignment solution over a square cost matrix
 * with integer valued costs.
 */
extern "C" SEXP ctc_integer_(
  const SEXP cost_  /*!< Pointer to cost matrix */
) {
  const Eigen::MatrixXi cost(
    Rcpp::as<Eigen::Map<Eigen::MatrixXi> >(cost_));
  std::vector<int> solution;
  abseil::min_cost_assignment::solve(cost, solution);
  Rcpp::IntegerVector x(solution.begin(), solution.end());
  return x;
};


/*
 * Find minimum assignment solution over a square cost matrix.
 * Values <= tol are taken equivalent.
 */
extern "C" SEXP ctc_double_(
  const SEXP cost_,  /*!< Pointer to cost matrix */
  const SEXP tol_    /*!< Pointer to scalar tolerance */
) {
  const Eigen::MatrixXd cost(
    Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(cost_));
  std::vector<int> solution;
  abseil::min_cost_assignment::tolerance( Rcpp::as<double>(tol_) );  
  abseil::min_cost_assignment::solve(cost, solution);
  Rcpp::IntegerVector x(solution.begin(), solution.end());
  return x;
};






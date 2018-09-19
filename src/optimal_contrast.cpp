// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// See http://www.w3.org/TR/WCAG20/#relativeluminancedef
inline double srgb_to_rgb(double channel) {
  double color = channel / 255;
  return color <= 0.03928 ? color / 12.92 : pow((color + 0.055) / 1.055, 2.4);
}

// See http://www.w3.org/TR/WCAG20/#getContrast-ratiodef
// [[Rcpp::export]]
double getContrast(NumericMatrix cand_rgb, NumericVector coloring, const arma::sp_mat& adj_mat, double beta) {
  NumericVector luminance(cand_rgb.rows());

  for (int i = 0; i< cand_rgb.rows(); i++) {
    luminance(i) =
      0.2126 * srgb_to_rgb(cand_rgb(i, 0)) +
      0.7152 * srgb_to_rgb(cand_rgb(i, 1)) +
      0.0722 * srgb_to_rgb(cand_rgb(i, 2));
  }

  double contrast = 0;
  int count = 0;

  for (arma::sp_mat::const_iterator i = adj_mat.begin(); i != adj_mat.end(); ++i) {
    double lum1 = luminance(coloring(i.row()) - 1);
    double lum2 = luminance(coloring(i.col()) - 1);

    double lum =
      lum1 > lum2 ?
        (lum1 + 0.05) / (lum2 + 0.05) :
        (lum2 + 0.05) / (lum1 + 0.05);

    contrast += pow(lum - 1 + 1e-3, -beta);
    count++;
  }

  return contrast / count;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
getContrast(matrix(1), 1, Matrix::sparseMatrix(1, 1, x = 1), 1.5)
*/

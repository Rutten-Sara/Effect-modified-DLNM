

Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> Hess_logpxi(const Eigen::SparseMatrix<double> &Qv,
                                        const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv) {
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;
  
  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }
  
  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  result += Qv;
  result *= -1;
  
  return result;
}
')





Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::SparseMatrix<double> Hess_logpxi_NB(const Eigen::SparseMatrix<double> &Qv,
                                          const Eigen::SparseMatrix<double> & Wnb,
                                        const Eigen::SparseMatrix<double> &MVnb,
                                        const Eigen::SparseMatrix<double> &Xv) {
                                     
   // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_M_Xv;
  result -= Xt_W_Xv;
  result -= Qv;
 
  return result;
 
}
')




Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::MatrixXd solve_sparse_cholesky(const Eigen::SparseMatrix<double> &Qv,
                                        const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv,
                                     const Eigen::VectorXd &b) {
                                     
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;
  
  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }
  
  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  result += Qv;
  result *= -1;
  
  result += 1e-5 * Eigen::MatrixXd::Identity(result.rows(), result.cols()).sparseView();


  // Step 4: Solve result*x = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(result);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");

  return solver.solve(b);
 
}
')




Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::MatrixXd solve_sparse_cholesky_NB(const Eigen::SparseMatrix<double> &Qv,
                                          const Eigen::SparseMatrix<double> & Wnb,
                                        const Eigen::SparseMatrix<double> &MVnb,
                                        const Eigen::SparseMatrix<double> &Xv,
                                     const Eigen::VectorXd &b) {
                                     
   // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_M_Xv;
  result -= Xt_W_Xv;
  result -= Qv;
  
  result += 1e-5 * Eigen::MatrixXd::Identity(result.rows(), result.cols()).sparseView();


  // Step 4: Solve result*x = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(result);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");

  return solver.solve(b);
 
}
')


Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double determinant_Pv(const Eigen::SparseMatrix<double> &Pv) {
     // Perform LDLT decomposition
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(Pv);

  if (ldlt.info() != Eigen::Success) {
    stop("LDLT decomposition failed");
  }

  // log-determinant = sum(log(diagonal of D))
  Eigen::VectorXd D_diag = ldlt.vectorD();  // for SimplicialLDLT, this is the diagonal of D
  if ((D_diag.array() <= 0).any()) {
    stop("Non-positive values in D; log-determinant is not defined");
  }

  double logdet = D_diag.array().log().sum();
  return logdet;

}
')





Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::SparseMatrix<double> XWX_func (const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv) {
                                     
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;
  
  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }
  
  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;
  
  // Result = Xt_W_Xv 
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  
  return result;
}
')


Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::SparseMatrix<double> XWX_func_NB (const Eigen::SparseMatrix<double> &Xv,
                  const Eigen::SparseMatrix<double> &Wnb) {
                                     
    // Construct matrix
    Eigen::SparseMatrix<double> Pv = Xv.transpose() * Wnb*Xv;

  
  return Pv;
}
')


# DIC calculation


calculate_DIC_tot <- function(xi_mode, Prec, Qv_mode, Xv, y, offset, v_mode, model.family = "poisson"){
  # DIC of mean
  
  if(model.family == "poisson"){
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      return(sum(y * log(Cvxi_mean) - Cvxi_mean-lfactorial(y)))
    }

  } else{
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      varval <-  Cvxi_mean + (1 / exp(v_mode[1])) * (Cvxi_mean ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_mean * (1/varval - (Cvxi_mean/(varval^2)) *
                                                 (1 + 2*Cvxi_mean*(1/exp(v_mode[1])))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_mean) ^ 2) * (1 / varval))
      gammanb <- exp(v_mode[1]) * log(Cvxi_mean / (Cvxi_mean+ exp(v_mode[1])))
      bgammanb <- (-1) * (exp(v_mode[1])^2) * log(exp(v_mode[1])/(exp(v_mode[1]) + Cvxi_mean))
      
      return( sum((1/exp(v_mode[1]))*((y * gammanb) - bgammanb) +
                         lgamma(y + exp(v_mode[1])) - lgamma(exp(v_mode[1])) - lgamma(y+1)))
    }

    
    
  }
  Dev_of_mean <- -2*log_lik(xi_mode, y, Xv, offset, v_mode)
  
  # pd
  pd_est <- dim(Prec)[1]-Hutchinson_pD(Prec, Qv_mode)
  
  DIC <- Dev_of_mean + 2*pd_est
  
  return(list(DIC = DIC, pd = pd_est))
}




Rcpp::sourceCpp(code = '
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
double Hutchinson_pD(
     const Eigen::SparseMatrix<double> &Ipost,
    const Eigen::SparseMatrix<double> &Ilike,
    int R = 5000)  // number of random probes
{
  const int p = Ipost.rows();
  double trace_est = 0.0;
  
  // Sparse Cholesky factorization of Ipost
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Ipost);
  if (chol.info() != Eigen::Success) stop("Cholesky decomposition failed");
  
  Rcpp::RNGScope scope;
  for (int r = 0; r < R; r++) {
    // ±1 Rademacher random vector
    Eigen::VectorXd z(p);
    for (int i=0; i<p; i++) z[i] = (R::runif(0,1) < 0.5 ? -1.0 : 1.0);

    // b = I_like * z
    Eigen::VectorXd b = Ilike * z;

    // Solve I_post * x = b
    Eigen::VectorXd x = chol.solve(b);

    trace_est += z.dot(x);
  }

  trace_est /= R;
  return trace_est;
}

')



library(irlba)
svd_k_irlba <- function(X, k){
  s <- irlba(X, nv = k, nu = k)
  list(U = s$u, S = Diagonal(x=s$d), V = s$v)
}


approx_XWX_from_svd <- function(U, S, V, w){
  # U: n x k, S: k x k, V: p x k, w: length n
  # compute Ut W U = t(U) %*% (w * U) efficiently
  Uw <- U * w            # scales rows of U by w (vector recycling)
  UtWU <- crossprod(U, Uw)   # k x k
  M <- S %*% UtWU %*% S      # k x k
  V %*% M %*% t(V)           # p x p approximate Gram
}




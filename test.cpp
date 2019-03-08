#include <RcppArmadillo.h>

using namespace Rcpp;
// using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat remove_col(arma::mat x, int which){
  x.shed_col(which);
  return(x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec select_col(arma::mat x, int which){
  arma::vec z = x.col(which);
  return(z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat select_row(arma::mat x, int which){
  arma::mat z = x.row(which);
  return(z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat select_row_remove_col(arma::mat x, int which_row, int which_col){
  arma::mat z = x.row(which_row);
  z.shed_col(which_col);
  return(z);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec offdiag_extract(arma::mat& A, int k) {
  return A.diag(k);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double select_row_select_col(arma::mat x, int which_row, int which_col){
  arma::vec z = x.row(which_row);
  double y = z[which_col];
  return(y);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat slice_extract(arma::cube X, int which){
   arma::mat temp = X.slice(which);
  return(temp);
  }


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat remove_row(arma::mat x, int which){
  x.shed_row(which);
  return(x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List test(arma::mat X, int samples, arma::mat int_mat, int temp) {
   
  // number of rows
   float n = X.n_rows;

   // number of columns
   int p = X.n_cols;
   
   // for intercept only design matrix
   arma::mat intercept_mat(n, 1, arma::fill::ones);
   
   // arma::mat mn(p, 1, arma::fill::ones);
   
   // matrix-F starting matrix
   // arma::mat Psi = arma::eye<arma::mat>(p, p);
   
   arma::mat b_inv = arma::eye<arma::mat>(p, p);
   
   arma::mat sigma_inv = arma::eye<arma::mat>(p, p);
   
   arma::mat check = arma::eye<arma::mat>(p, p);
   
   arma::mat vb = arma::eye<arma::mat>(p, p);
   // residual correlation matrix

  
   // ordinal levels
   arma::vec n_levels = unique(X.col(0));
   
  Rcpp::Range thresh_sampled =  seq(min(X.col(0)) + 1, max(X.col(0)) - 1 );
   
   // int i;  int j;
   arma::mat mattemp = arma::mat(n,p);
 
   arma::mat e(p, 1, arma::fill::zeros);

    arma::mat mattemp2(n,p);
   arma::colvec mattemp3(p);
   
   arma::mat mat2(p,p);
   // arma::mat  vt(n, p);
   
 
   arma::colvec mattemp4(p);
   
   
   arma::rowvec b2(p, arma::fill::zeros); 
   
  // note 0 = 1; p - 1 give 0:(p-1)--i.e., 1:p (in R)
  arma::mat collect(samples, 10);
  
  arma::mat collect_R(samples, 10);
  
  // arma::mat  mm(p,1);
  
 
  
  arma::cube  Psi = arma::cube(p, p, 1, arma::fill::zeros);
  arma::cube cor_save = arma::cube(p,p, samples);
  
 
  
  
   
 
  
  arma::mat R_top(n,p);
   
   arma::mat  R_bottom = arma::mat(n,p);
   
   // arma::rowvec mattemp4(p);
   
  
  for(int i = 0; i < p; i ++){
    Psi(i,i,0) = 1;
}
 
 
 
 
 
 
  arma::mat slice_thresh = arma::mat(n,p);
  
  arma::mat  b(p, 1, arma::fill::zeros);
 
  
   arma::mat ss = arma::mat(p, 1);
  
 
   arma::cube  thresh_mat = arma::zeros<arma::cube>(p, n_levels.size() + 1, samples);
   
   arma::cube c_thresh_mat = arma::cube(p, n_levels.size() + 1, samples, arma::fill::zeros); 
   
   
   arma::cube cor_mat = arma::cube(p,p,1, arma::fill::zeros);
   
   for(int i = 0; i < p; i ++){
     for(int j = 0; j < p; j ++){
       
       cor_mat(i,j,0) = int_mat(i,j);
     }
   }
   
   
   
   // arma::cube zstar(n, p, samples, arma::fill::zeros); 
   arma::mat zstar(n, p, arma::fill::zeros); 
   
   int arr[] = {0, 5};

for(int k = 0; k < samples; k ++){
   
   for(int i = 0; i < p; i ++ ){
     for(int j = 0; j < 2; j ++){
       
       if(j == 0){

             c_thresh_mat(i,arr[j],k) = -arma::datum::inf;

             thresh_mat(i,arr[j], k) =  -arma::datum::inf;
       } if(j == 1){

         c_thresh_mat(i,arr[j],k) = arma::datum::inf;

         thresh_mat(i, arr[j],k) = arma::datum::inf;
       }
     }
   }

}

  for(int k = 0; k < samples; k ++){

    arma::mat m = intercept_mat * trans(b.col(0));
    
    for(int i = 0; i < max(X.col(0)) - 1; i ++){

      for(int j = 0; j <= p - 1; j ++){

              thresh_mat(j, i + 2, k) = R::qnorm5(sum(X.col(j) <= thresh_sampled[i] ) / n,
                                       
                                      - R::qnorm5(sum(X.col(j) == 1) / n, 0, 1, TRUE, FALSE), 1, TRUE, FALSE);

      }
  }


   arma::mat  mm(n,p);



   for(int i = 0; i < p ; i++){


   // row i and remove column i--i.e., mat[ i, -i]
   arma::mat mat_i_not_i = select_row_remove_col(cor_mat.slice(0), i, i);


   // remove row i and remove column i--i.e., mat[ - i, -i]
   arma::mat mat_not_i_not_i = inv(remove_row(remove_col(cor_mat.slice(0), i), i));


   // updated: subset of multivariate location
   arma::colvec mm = trans(trans(m.col(i)) + mat_i_not_i * mat_not_i_not_i *  trans(remove_col(zstar, i) - remove_col(m, i)));


   // updated: subset of multivariate scale
   arma::mat ss =  select_row(cor_mat.slice(0).col(i) , i) -  mat_i_not_i  * mat_not_i_not_i * trans(mat_i_not_i);


  for(int j = 0; j < n ; j ++ ){
      zstar(j, i) = R::qnorm(R::runif(
                                     // minimum
                                     R::pnorm(thresh_mat.slice(0)(i , (X.col(i)[j] - 1)), mm(j), sqrt(ss(0)), TRUE, FALSE),
                                    // maximum
                                     R::pnorm(thresh_mat.slice(0)(i , X.col(i)[j]), mm(j), sqrt(ss(0)), TRUE, FALSE)),
                                    // location and scale
                                     mm(j), sqrt(ss(0)), TRUE, FALSE) ;


    }
  }



  
     
 for(int i = 0; i < max(X.col(0)) - 1; i ++){

   // note 0 = 1; 1 - max gives largest threshold
   for(int j = 0; j < p ; j ++){

     arma::mat slice_c_thresh = thresh_mat.slice(0);

     c_thresh_mat(j, thresh_sampled[i], 0) = R::qnorm5(R::runif(R::pnorm(0,  slice_c_thresh(j , thresh_sampled[i]), .01, TRUE, FALSE), 1),
                                                       slice_c_thresh(j , thresh_sampled[i]),   0.01, TRUE, FALSE);
       // location and scale, TRUE, FALSE);

   }
 }

 


for(int i = 0; i < p ; i ++){
  for(int j = 0; j < n ; j ++ ){
   
      R_top(j, i) =  (R::pnorm5(c_thresh_mat(i , X.col(i)[j], 0) - mm(j), 0, 1, TRUE, FALSE) -
                      R::pnorm5(c_thresh_mat(i , X.col(i)[j] - 1, 0) - mm(j), 0, 1, TRUE, FALSE));
                                
     R_bottom(j, i) = (R::pnorm5(thresh_mat(i , X.col(i)[j], 0) - mm(j), 0, 1, TRUE, FALSE) -
                       R::pnorm5(thresh_mat(i , X.col(i)[j] - 1, 0) - mm(j), 0, 1, TRUE, FALSE));
                                                              
                                                              
  }
}

   for(int i = 0; i < p; i ++){
    mattemp4(i) = accu(log(R_top.col(i)) - log(R_bottom.col(i))) ;

   }


           for(int j = 0; j < p; j ++){
           if( mattemp4[j] > log(R::runif(0, 1) )){

              thresh_mat.slice(0).row(j) = c_thresh_mat.slice(0).row(j);
       
     }  
   
  }

             arma::mat vb(inv(inv(cor_mat.slice(0) / n)));
             
             
             
             arma::colvec mn = arma::colvec(vb * trans((intercept_mat.t() * zstar *  trans(inv(cor_mat.slice(0))))));
          
             arma::mat b_temp = mn.t() + (b.randn().t() * chol(vb));
             
             for(int i = 0; i < p; i ++){
               b(i) = b_temp(i) ; 
               
             }
             

             arma::mat e = zstar - (trans(b_temp.t() * intercept_mat.t()));
             arma::mat  vt =  trans(e) * e ;
             
             
            arma::mat  sig_inv = wishrnd(inv(vt + Psi.slice(0)), n + 10 - 1);
            
           
             
             arma::vec temp6 =  1 / sqrt(sig_inv.diag());
        
             cor_save.slice(k) = diagmat(temp6) * sig_inv * diagmat(temp6) * - 1;
  
             arma::mat cov_mat = inv(sig_inv);
             
             arma::vec temp5 =  1 / sqrt(cov_mat.diag());
             
             cor_mat.slice(0) = diagmat(temp5) * cov_mat * diagmat(temp5);
             
             Psi.slice(0) = wishrnd(inv(sig_inv + b_inv * 1000), 1000 + p - 2);
  
  
  }
 

 List ret;
 ret["thresh_samp"] = thresh_sampled;
 ret["zstar"] = zstar;
 ret["c_thresh"] = c_thresh_mat;
 ret["thresh"] = thresh_mat;
 ret["n"] = n;
 ret["Psi"] = Psi;
 ret["R_top"] = R_top;
 ret["R_bottom"] = R_bottom;
 ret["mattemp4"] = mattemp4;
 ret["b"] = b;
 ret["cor_mat"] = cor_mat;
 ret["cor_save"] = cor_save;
 
 //d

 return ret;
 
}
 



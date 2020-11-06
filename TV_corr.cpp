#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// Functions


// Create a block-diagonal matrix
mat blockDiag(const mat& X, const mat& Y){
  mat A = zeros(X.n_rows, Y.n_cols);
  mat B = zeros(Y.n_rows, X.n_cols);
  mat Z = join_cols(join_rows(X,A),join_rows(B,Y));
  return Z;
} 


// Inverse of a matrix based on the SVD decomposition
mat invSVD(const mat& X) {
  mat U, V;
  vec S;
  svd(U, S, V, X, "standard");
  return V*diagmat(1/S)*U.t();
}


// Link function that bounds its argument between -1 and 1
double Clink(double& x){
  double y;
  y = x/sqrt(1+pow(x,2));
  return y;
}


// Append a scalar to a vector
vec vecscal(const vec& x, double& y){
  vec z(x.size()+1);
  z.subvec(0,x.size()-1) = x;
  z(x.size()) = y;
  return z;
}


// Convert the element of a list into a matrix
mat list2mat(Rcpp::List& y, int& j){
  mat x = y[j];
  return x;
}


// Convert the element of a list into a vector
vec list2vec(Rcpp::List& y, int& j){
  vec x = y[j];
  return x;
}


// Stratified resampling
vec stratifiedResampling_rcpp(vec& w, const int& M){
  w = w/sum(w);
  vec cdf_w = cumsum(w);
  vec resample_set = zeros(M);
  double u;
  for (int m = 0; m < M; m++){
    u = m/double(M) + ((m+1)/double(M)-m/double(M))*as_scalar(randu(1));
    if (any(cdf_w < u)){
      resample_set(m) = index_max(find(cdf_w < u)) + 1; // length(cdf_w[which(cdf_w < runif(1,(m-1)/M,m/M))])+1;
    } else {
      resample_set(m) = 0;
    }
  }
  return resample_set;
}


// Extract the columns of a matrix according to a vector of column-indeces, and store them into a new matrix
mat repcol(const mat& X, const vec& ind) {
  mat Xsub = X;
  int l = ind.size();
  for (int j = 0; j < l; j++){
    Xsub.col(j) = X.col(ind(j));
  }
  return Xsub;
}


// Extract the elements of a list according to a vector of element-indeces, and store them into a new list
Rcpp::List replist(const Rcpp::List& X, const vec& ind) {
  Rcpp::List Xsub = X;
  int l = ind.size();
  for (int j = 0; j < l; j++){
    Xsub[j] = X[ind(j)];
  }
  return Xsub;
}


// Column means of a matrix
vec colMeans(const mat& A){
  int l = A.n_cols;
  vec b = zeros(l);
  for (int j = 0; j < l; j++){
    b(j) = mean(A.col(j));
  } 
  return b;
}


// Kalman filter estimation of the cubic splines model
double KF_CC_splines_rcpp(const vec& par, const mat& y, const mat& se, const bool& opti, const bool& outofsample, const double& parP10, 
                          const int& nstates, int& d, const int& k, const mat& W, const bool& hyper_tan, const bool& restricted){  
  int len = y.n_cols;
  double delta = par(11);
  mat sd_nu = diagmat(exp(par.subvec(3,7))); 
  colvec x10 = zeros(nstates);
  Rcpp::List Pttm1(len+1);
  Pttm1[0] = blockDiag(parP10*diagmat(ones(17)), diagmat(ones(1)));
  Pttm1[0] = blockDiag(Pttm1[0], (1-pow(delta,2))*diagmat(ones(4)));
  Pttm1[0] = blockDiag(Pttm1[0], diagmat(ones(1)));
  Pttm1[0] = blockDiag(Pttm1[0], (1-pow(delta,2))*diagmat(ones(3)));
  Pttm1[0] = blockDiag(Pttm1[0], diagmat(ones(1)));
  Pttm1[0] = blockDiag(Pttm1[0], (1-pow(delta,2))*diagmat(ones(3)));
  Pttm1[0] = blockDiag(Pttm1[0], parP10*diagmat(ones(nstates-30)));
  Rcpp::List Ptt(len);
  Ptt[0] = zeros(nstates,nstates); 
  mat xtt = zeros(nstates,len); 
  mat xttm1 = zeros(nstates,len+1); 
  xttm1.col(0) = x10;
  colvec epshatoutofsample = zeros(y.n_rows);
  colvec epshatinsample = zeros(y.n_rows);
  mat Fmatrix = zeros(y.n_rows,y.n_rows);
  mat Fmatrixinv = zeros(y.n_rows,y.n_rows);
  mat H = blockDiag(zeros(5,5), diagmat(exp(2*par(10))*ones(1)));
  mat Kg = zeros(nstates,y.n_rows);
  mat Zy;
  mat Zx;
  mat Z;  
  d = d-1;  
  mat R = diagmat(ones(nstates));
  mat Q;
  mat D = blockDiag(diagmat(zeros(1)), exp(par(0))*diagmat(ones(1)));
  D = blockDiag(D, exp(par(1))*diagmat(ones(11)));
  D = blockDiag(D, exp(par(2))*diagmat(ones(4)));
  D = blockDiag(D, sd_nu);
  D = blockDiag(D, diagmat(zeros(9)));
  D = blockDiag(D, exp(par(8))*diagmat(ones(1)));
  D = blockDiag(D, exp(par(9))*diagmat(ones(11)));
  double rho;
  double gamma;  
  double pi = 3.141592653589793238462643383280;

  mat Tymu = diagmat(ones(2));
  Tymu(0,1) = 1;
  Rcpp::List C(5);
  mat Cmat;
  for (int l = 0; l < 5; l++){
    Cmat = diagmat(cos((pi*(l+1))/6)*ones(2));
    Cmat(0,1) = sin((pi*(l+1))/6);
    Cmat(1,0) = -sin((pi*(l+1))/6);
    C[l] = Cmat;
  }
  mat Tyomega = blockDiag(C[0],C[1]);
  Tyomega = blockDiag(Tyomega,C[2]);
  Tyomega = blockDiag(Tyomega,C[3]);
  Tyomega = blockDiag(Tyomega,C[4]);
  Tyomega = blockDiag(Tyomega,-1*diagmat(ones(1)));
  mat Tylambda = diagmat(ones(4));
  mat TyE = join_rows(join_cols(zeros(9,4), diagmat(ones(4))), zeros(13,1));
  TyE = join_rows(TyE, join_cols(zeros(1,4), delta*diagmat(ones(4)), zeros(8,4)));
  TyE = join_rows(TyE, join_cols(zeros(5,4), diagmat(ones(4)), zeros(4,4)));
  mat Ty = blockDiag(Tymu, Tyomega);
  Ty = blockDiag(Ty, Tylambda);
  Ty = blockDiag(Ty, TyE);
  mat Tx = blockDiag(Tymu, Tyomega);
  mat Tmatrix = blockDiag(Ty, Tx);

  double logl = 0;

  for (int i = 0; i < len; i++){
    Zy = ones(1,2);
    Zy(0,1) = 0;
    Zy = repmat(Zy, 1, 6);
    Zy = join_rows(Zy, diagmat(ones(1)));
    Zy = join_cols(Zy, Zy);
    Zy = join_cols(Zy, Zy);
    Zy = join_cols(Zy,Zy.submat(0,0,0,Zy.n_cols-1));
    Zy = join_rows(Zy, join_cols(zeros(1,4), diagmat(ones(4))));
    Zy = join_rows(Zy, join_rows(diagmat(se.row(i)), zeros(5,8)));
    Zx = ones(1,2);
    Zx(0,1) = 0;
    Zx = repmat(Zx, 1, 6);
    Zx = join_rows(Zx, diagmat(ones(1)));
    Z = blockDiag(Zy,Zx);  
    
    epshatoutofsample = y.col(i) - Z*xttm1.col(i);
    Fmatrix = Z*list2mat(Pttm1,i)*Z.t() + H;
    Fmatrixinv = invSVD(Fmatrix);
    Kg = Tmatrix*list2mat(Pttm1,i)*Z.t()*Fmatrixinv;
    xtt.col(i) = xttm1.col(i)+list2mat(Pttm1,i)*Z.t()*Fmatrixinv*epshatoutofsample;
    epshatinsample = y.col(i)-Z*xtt.col(i);
    Ptt[i] = list2mat(Pttm1,i)-list2mat(Pttm1,i)*Z.t()*Fmatrixinv*Z*list2mat(Pttm1,i);
    
    if (restricted){
      gamma = par(par.size()-k);
    } else {
      gamma = as_scalar(W.submat(i,0,i,(k-1))*par.subvec(par.size()-k,par.size()-1));
    }
    
    if (hyper_tan){
      rho = tanh(gamma);
    } else {
      rho = Clink(gamma);
    }
    
    R(31,1) = rho;
    R(1,31) = rho;
    Q = D*R*D;
    Pttm1(i+1) = Tmatrix*list2mat(Pttm1,i)*(Tmatrix-Kg*Z).t()+Q;
    xttm1.col(i+1) = Tmatrix*xttm1.col(i) + Kg*epshatoutofsample;
    
    if (outofsample) {
      if (i <= d ){
        logl = logl - (y.n_rows*log(2*pi))/2;      
      } else {
        logl = logl - (y.n_rows*log(2*pi))/2 - log(det(Fmatrix))/2 - as_scalar(epshatoutofsample.t()*Fmatrixinv*epshatoutofsample)/2;
      }
    } else {
      if (i <= d ){
        logl = logl - (y.n_rows*log(2*pi))/2;
      } else {
        logl = logl - (y.n_rows*log(2*pi))/2 - log(det(Fmatrix))/2 - as_scalar(epshatinsample.t()*Fmatrixinv*epshatinsample)/2;
      }
    }
  }
  
  logl = -logl;
  
  return logl;
}  


// [[Rcpp::export]]
Rcpp::List ucminf_rcpp_splines(const vec& init_val, const mat& y, const mat& se, const bool& opti, const bool& outofsample, const double& parP10, 
                               const int& nstates, int& d, const int& k, const mat& W, const bool& restricted, const bool& hyper_tan, Rcpp::List& control){ 
  // Extract R's ucminf function
  Rcpp::Environment stats("package:ucminf"); 
  Rcpp::Function ucminf = stats["ucminf"];
  
  // Call the ucminf function from R in C++ 
  Rcpp::List opt_results = ucminf(Rcpp::_["par"]    = init_val,
                                  // Make sure this function is not exported!
                                  Rcpp::_["fn"]     = Rcpp::InternalFunction(&KF_CC_splines_rcpp),
                                  Rcpp::_["y"] = y,
                                  Rcpp::_["se"] = se,
                                  Rcpp::_["opti"] = opti,
                                  Rcpp::_["outofsample"] = outofsample,
                                  Rcpp::_["parP10"] = parP10,
                                  Rcpp::_["nstates"] = nstates,
                                  Rcpp::_["d"] = d,
                                  Rcpp::_["k"] = k,
                                  Rcpp::_["W"] = W,
                                  Rcpp::_["hyper_tan"] = hyper_tan,
                                  Rcpp::_["restricted"] = restricted,
                                  Rcpp::_["control"] = control);
  
  // Extract out the estimated parameter values
  vec est_par = Rcpp::as<arma::vec>(opt_results[0]);      // estimated parameters
  double value = opt_results[1];      // value of the log-likelihood at the estimates
  
  // Return estimated parameters and log-likelihood value
  Rcpp::List ret;
  ret["par"] = est_par;
  ret["value"] = value;
  return ret;
}


// Compute the prediction step of the Kalman filter
// [[Rcpp::export]]
Rcpp::List KF_t_rcpp(const vec& par, const vec& y, const vec& se, vec& xttm1, mat& Pttm1, const bool& hyper_tan, double& gamma_draw, const int& nstates){  
  vec epshatoutofsample;
  mat Fmatrix;
  mat Fmatrixinv;
  mat H = blockDiag(zeros(5,5), diagmat(exp(2*par(10))*ones(1)));
  mat Kg;
  vec xtt;
  vec epshatinsample;
  mat Ptt;
  mat sd_nu = diagmat(exp(par.subvec(3,7)));
  mat R = diagmat(ones(nstates));
  mat Q;
  mat D = blockDiag(diagmat(zeros(1)), exp(par(0))*diagmat(ones(1)));
  D = blockDiag(D, exp(par(1))*diagmat(ones(11)));
  D = blockDiag(D, exp(par(2))*diagmat(ones(4)));
  D = blockDiag(D, sd_nu);
  D = blockDiag(D, diagmat(zeros(9)));
  D = blockDiag(D, exp(par(8))*diagmat(ones(1)));
  D = blockDiag(D, exp(par(9))*diagmat(ones(11)));
  double rho; 
  double delta = par(11);   
  double logl;  
  double pi = 3.141592653589793238462643383280;

  mat Tymu = diagmat(ones(2));
  Tymu(0,1) = 1;
  Rcpp::List C(5);
  mat Cmat;
  for (int l = 0; l < 5; l++){
    Cmat = diagmat(cos((pi*(l+1))/6)*ones(2));
    Cmat(0,1) = sin((pi*(l+1))/6);
    Cmat(1,0) = -sin((pi*(l+1))/6);
    C[l] = Cmat;
  }
  mat Tyomega = blockDiag(C[0],C[1]);
  Tyomega = blockDiag(Tyomega,C[2]);
  Tyomega = blockDiag(Tyomega,C[3]);
  Tyomega = blockDiag(Tyomega,C[4]);
  Tyomega = blockDiag(Tyomega,-1*diagmat(ones(1)));
  mat Tylambda = diagmat(ones(4));
  mat TyE = join_rows(join_cols(zeros(9,4), diagmat(ones(4))), zeros(13,1));
  TyE = join_rows(TyE, join_cols(zeros(1,4), delta*diagmat(ones(4)), zeros(8,4)));
  TyE = join_rows(TyE, join_cols(zeros(5,4), diagmat(ones(4)), zeros(4,4)));
  mat Ty = blockDiag(Tymu, Tyomega);
  Ty = blockDiag(Ty, Tylambda);
  Ty = blockDiag(Ty, TyE);
  mat Tx = blockDiag(Tymu, Tyomega);
  mat Tmatrix = blockDiag(Ty, Tx);

  mat Zy;
  mat Zx;
  mat Z;
  Zy = ones(1,2);
  Zy(0,1) = 0;
  Zy = repmat(Zy, 1, 6);
  Zy = join_rows(Zy, diagmat(ones(1)));
  Zy = join_cols(Zy, Zy);
  Zy = join_cols(Zy, Zy);
  Zy = join_cols(Zy,Zy.submat(0,0,0,Zy.n_cols-1));
  Zy = join_rows(Zy, join_cols(zeros(1,4), diagmat(ones(4))));
  Zy = join_rows(Zy, join_rows(diagmat(se), zeros(5,8)));
  Zx = ones(1,2);
  Zx(0,1) = 0;
  Zx = repmat(Zx, 1, 6);
  Zx = join_rows(Zx, diagmat(ones(1)));
  Z = blockDiag(Zy,Zx);

  epshatoutofsample = y - Z*xttm1;
  Fmatrix = Z*Pttm1*Z.t() + H;
  Fmatrixinv = invSVD(Fmatrix);
  Kg = Tmatrix*Pttm1*Z.t()*Fmatrixinv;
  xtt = xttm1+Pttm1*Z.t()*Fmatrixinv*epshatoutofsample;
  epshatinsample = y-Z*xtt;
  Ptt = Pttm1-Pttm1*Z.t()*Fmatrixinv*Z*Pttm1;
  
  if (hyper_tan){
    rho = tanh(gamma_draw);
  } else {
    rho = Clink(gamma_draw);
  }
  
  R(31,1) = rho;
  R(1,31) = rho;
  Q = D*R*D;
  
  Pttm1 = Tmatrix*Pttm1*(Tmatrix-Kg*Z).t()+Q;
  xttm1 = Tmatrix*xttm1 + Kg*epshatoutofsample;
  
  logl = - (y.size()*log(2*pi))/2 - log(det(Fmatrix))/2 - as_scalar(epshatoutofsample.t()*Fmatrixinv*epshatoutofsample)/2;
  
  Rcpp::List ret;
  ret["value"] = -logl;
  ret["xtt"] = xtt;
  ret["xttm1"] = xttm1;
  ret["Pttm1"] = Pttm1;
  ret["Ptt"] = Ptt;
  return ret;
}


// Rao-Blackwellised bootstrap filter estimation of the state vector of the nonlinear model
// [[Rcpp::export]]
Rcpp::List boot_filter_CC_rcpp(const int& draw_m, const vec& tau_hat, const mat& y, const mat& se, const int& nstates, const bool& hyper_tan,
                               const mat& Rsel, const vec& states_noerr, const double& init_gamma){
  int len = y.n_cols;
  Rcpp::List KF_results;
  mat alpha_t = zeros(nstates+1, draw_m);
  vec w_tilde = zeros(draw_m);
  vec p_density = zeros(draw_m);
  double p_den_m;
  vec ESS = zeros(len);
  vec CV = zeros(len);
  vec w_t = zeros(draw_m);
  vec w_t_power = zeros(draw_m);

  mat Zymatrix = ones(1,1);
  mat Zxmatrix = Zymatrix;
  mat Zmatrix = blockDiag(Zymatrix,Zxmatrix);

  vec att_BF = zeros(len);
  vec Ptt_BF = zeros(len);
  Rcpp::List Pttm1_KF(draw_m);
  Rcpp::List Pttm1_KF_next(draw_m);
  Rcpp::List Ptt_KF(len);
  mat att_KF = zeros(nstates, len);

  double gamma_t;
  double rho_t;
  vec level_t;
  vec alpha_tm;
  mat Pttm1_tm;
  double gamma_tm;
  vec resample_set;

  mat R = diagmat(ones(nstates-states_noerr.size()));
  mat Q;
  mat sd_nu = diagmat(exp(tau_hat.subvec(3,7)));
  mat D = blockDiag(exp(tau_hat(0))*diagmat(ones(1)), exp(tau_hat(1))*diagmat(ones(11)));
  D = blockDiag(D, exp(tau_hat(2))*diagmat(ones(4)));
  D = blockDiag(D, sd_nu);
  D = blockDiag(D, exp(tau_hat(8))*diagmat(ones(1)));
  D = blockDiag(D, exp(tau_hat(9))*diagmat(ones(11)));
  mat Dchol = chol(D);
  mat M;

  // Initialisation
  for (int m = 0; m < draw_m; m++){
    gamma_t = as_scalar(exp(tau_hat(tau_hat.size()-1))*randn(1) + init_gamma);      // initialise gamma_t

    if (hyper_tan){
      rho_t = tanh(gamma_t);
    } else {
      rho_t = Clink(gamma_t);
    }

    R(21,0) = rho_t;
    R(0,21) = rho_t;
    Q = D*R*D;
    M = chol(Q);
    level_t = Rsel*M.t()*randn(nstates-states_noerr.size(),1);
    alpha_t.col(m) = vecscal(level_t, gamma_t);
    w_t = ones(draw_m)/draw_m;
    Pttm1_KF[m] = Rsel*Q*Rsel.t();
  }

  // Loop over time
  for (int i = 0; i < len; i++){
    for (int m = 0; m < draw_m; m++){

      alpha_tm = vectorise(alpha_t.submat(0,m,alpha_t.n_rows-2,m));
      Pttm1_tm = list2mat(Pttm1_KF,m);
      gamma_tm = as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m));

      KF_results = KF_t_rcpp(tau_hat.subvec(0,tau_hat.size()-2),vectorise(y.col(i)),vectorise(se.row(i)),alpha_tm,Pttm1_tm,hyper_tan,gamma_tm,nstates);

      // Compute weights
      p_den_m = KF_results["value"];
      p_density(m) = exp(-p_den_m);
      if (p_density(m) == 0){
        p_density(m) = 1e-20;
      }
      w_tilde(m) = w_t(m)*p_density(m);
    }

    // Normalised weights
    w_t = w_tilde/sum(w_tilde);

    ESS(i) = as_scalar(1/(sum(pow(w_t,2))));
    CV(i) = sqrt(mean(pow(draw_m*w_t-1,2)));

    // Resample
    w_t_power = pow(w_t,(log(sqrt(draw_m-1)/CV(i))));
    resample_set = stratifiedResampling_rcpp(w_t_power, draw_m);
    alpha_t = repcol(alpha_t,resample_set);
    Pttm1_KF = replist(Pttm1_KF,resample_set);

    w_t = ones(draw_m)/draw_m;

    // Estimate state variable and variance
    for (int m = 0; m < draw_m; m++){
      att_BF(i) = att_BF(i) + w_t(m)*as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m));
      Ptt_BF(i) = Ptt_BF(i) + w_t(m)*pow(as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m)),2);
    }
    Ptt_BF(i) = Ptt_BF(i) - pow(att_BF(i),2);

    if (i < len-1){
      // Regenerate data
      for (int m = 0; m < draw_m; m++){
        gamma_t = as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m)) + as_scalar(exp(tau_hat(tau_hat.size()-1))*randn(1));

        alpha_tm = vectorise(alpha_t.submat(0,m,alpha_t.n_rows-2,m));
        Pttm1_tm = list2mat(Pttm1_KF,m);

        KF_results = KF_t_rcpp(tau_hat.subvec(0,tau_hat.size()-2),vectorise(y.col(i)),vectorise(se.row(i)),alpha_tm,Pttm1_tm,hyper_tan,gamma_t,nstates);

        level_t = Rcpp::as<arma::vec>(KF_results["xttm1"]);
        alpha_t.col(m) = vecscal(level_t, gamma_t);
        Pttm1_KF[m] = KF_results["Pttm1"];

      }
    }

  }

  Rcpp::List ret;
  ret["att_BF"] = att_BF;
  ret["Ptt_BF"] = Ptt_BF;
  ret["p_density"] = p_density;
  ret["CV"] =  CV;
  ret["ESS"] = ESS;
  ret["resample_set"] = resample_set;
  return ret;
}

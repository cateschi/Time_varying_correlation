#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// the order in which you put the elements in the input matters!!
// is a and b are integers, the divion works as a/double(b)


mat blockDiag(const mat& X, const mat& Y){
  mat A = zeros(X.n_rows, Y.n_cols);
  mat B = zeros(Y.n_rows, X.n_cols);
  mat Z = join_cols(join_rows(X,A),join_rows(B,Y));
  return Z;
} 


mat invSVD(const mat& X) {
  mat U, V;
  vec S;
  svd(U, S, V, X, "standard");
  return V*diagmat(1/S)*U.t();
}


double Clink(double& x){
  double y;
  y = x/sqrt(1+pow(x,2));
  return y;
}


vec vecscal(const vec& x, double& y){
  vec z(x.size()+1);
  z.subvec(0,x.size()-1) = x;
  z(x.size()) = y;
  return z;
}


mat list2mat(Rcpp::List& y, int& j){
  mat x = y[j];
  return x;
}


vec list2vec(Rcpp::List& y, int& j){
  vec x = y[j];
  return x;
}


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


mat repcol(const mat& X, const vec& ind) {
  mat Xsub = X;
  int l = ind.size();
  for (int j = 0; j < l; j++){
    Xsub.col(j) = X.col(ind(j));
  }
  return Xsub;
}


Rcpp::List replist(const Rcpp::List& X, const vec& ind) {
  Rcpp::List Xsub = X;
  int l = ind.size();
  for (int j = 0; j < l; j++){
    Xsub[j] = X[ind(j)];
  }
  return Xsub;
}


vec colMeans(const mat& A){
  int l = A.n_cols;
  vec b = zeros(l);
  for (int j = 0; j < l; j++){
    b(j) = mean(A.col(j));
  } 
  return b;
}




double KF_CC_splines_rcpp(const vec& par, const mat& y, const mat& se, const bool& opti, const bool& outofsample, const double& parP10, 
                          const int& nstates, int& d, const int& k, const mat& W, const bool& hyper_tan, const bool& restricted){
  
  int len = y.n_cols;
  double delta = par(11);
  mat sd_nu = diagmat(exp(par.subvec(3,7)));
  //double sigma_Ry = par(0);
  //double sigma_Rx = par(1);
  //double sd_nu_y = par(2);
  //double sd_nu_x = par(3);
  
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
  
  //Build T;
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
  
  //initialization of loglikelihood
  double logl = 0;
  
  //KF recursions
  for (int i = 0; i < len; i++){
    
    //Bulild Z:
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
      //gamma = as_scalar(par(par.size()-k) + W.submat(i,1,i,(k-1))*par.subvec(par.size()-k+1,par.size()-1));
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
    
    
    //The optimization criterion
    if (outofsample) {
      if (i <= d ){
        logl = logl - (y.n_rows*log(2*pi))/2;      // in rcpp you cannot pre-multiply by 1/2*
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
  // Return log-likelihood
  return logl;
}  



// [[Rcpp::export]]
Rcpp::List ucminf_rcpp_splines(const vec& init_val, const mat& y, const mat& se, const bool& opti, const bool& outofsample, const double& parP10, 
                               const int& nstates, int& d, const int& k, const mat& W, const bool& restricted, const bool& hyper_tan, Rcpp::List& control){
  
  // Extract R's optim function
  Rcpp::Environment stats("package:ucminf"); 
  Rcpp::Function ucminf = stats["ucminf"];
  
  // Call the optim function from R in C++ 
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
  
  // Return estimated values
  Rcpp::List ret;
  ret["par"] = est_par;
  ret["value"] = value;
  return ret;
}




double IndInf_CC_rcpp(const vec& par, const int& sim_h, const vec& beta_hat, const mat& eps_sim, double& VARgamma, 
                      const bool& opti, const bool& outofsample, const double& parP10, int& d, const mat& W,
                      const bool& restricted, const bool& hyper_tan, const mat& se, const vec& states_noerr, const vec& init_val_CC,
                      const mat& Rsel, const vec& Msel, const vec& init_valII, Rcpp::List& control){
  
  int nstates = 43;
  int nvar = 6;
  int len = eps_sim.n_cols/double(sim_h);
  int k = W.n_cols;
  
  //double sigma_Ry = par(0);
  //double sigma_omegay = par(1);
  //double sigma_lambda = par(2);
  mat sd_nu = diagmat(exp(par.subvec(3,7)));
  //double sigma_Rx = par(8);
  //double sigma_omegax = par(9);
  mat eta;
  mat Rmat = diagmat(ones(nstates-states_noerr.size()));
  
  mat betatilde = zeros(sim_h, par.size());
  vec gammaII;
  vec rhoII(len);
  Rcpp::List Q(len);
  Rcpp::List R(len);
  Rcpp::List M(len);
  mat D = blockDiag(exp(par(0))*diagmat(ones(1)), exp(par(1))*diagmat(ones(11)));
  D = blockDiag(D, exp(par(2))*diagmat(ones(4)));
  D = blockDiag(D, sd_nu);
  D = blockDiag(D, exp(par(8))*diagmat(ones(1)));
  D = blockDiag(D, exp(par(9))*diagmat(ones(11)));
  mat Dchol = chol(D);
  
  mat alpha = zeros(nstates,len);
  vec nu;
  mat y = zeros(nvar,len);
  mat Zy;
  mat Zx;
  mat Z;
  
  double delta = par(11);
  
  double pi = 3.141592653589793238462643383280;
  
  //Build T;
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
  
  Rcpp::List objoptII(2);
  vec estpar;
  double var_est_gamma;
  double opt_fun;
  vec crit_min;
  
  for (int h = 0; h < sim_h; h++){
    
    //generate TV corr according to a random walk model:
    gammaII = cumsum(exp(par(par.size()-1))*vectorise(eps_sim.submat(1,(h+1)*len-len,1,(h+1)*len-1)));
    if (hyper_tan){
      for (int i = 0; i < len; i++){
        rhoII(i) = tanh(gammaII(i));
      }
    } else {
      for (int i = 0; i < len; i++){
        rhoII(i) = Clink(gammaII(i));
      }
    }
    
    //generate Innovations of states:
    eta = eps_sim.submat(2,(h+1)*len-len,eps_sim.n_rows-1,(h+1)*len-1);
    for (int i = 0; i < len; i++){
      Rmat(0,21) = rhoII(i);
      Rmat(21,0) = rhoII(i);
      R[i] = Rmat;
      Q[i] = D*Rmat*D;
      M[i] = Dchol*chol(Rmat)*Dchol;
      eta.col(i) = list2mat(M,i).t()*eta.col(i);
    }
    
    //generate state variables
    alpha.col(0) = Rsel*eta.col(0);
    for (int i = 1; i < len; i++){
      alpha.col(i) = Tmatrix*alpha.col(i-1) + Rsel*eta.col(i);
    }
    
    //Generate the innovations in the observation equation
    nu = exp(par(10))*vectorise(eps_sim.submat(0,(h+1)*len-len,0,(h+1)*len-1));
    
    
    //Generate the observable variables - the final model
    for (int i = 0; i < len; i++){
      
      //Bulild Z:
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
      
      y.col(i) = Z*alpha.col(i) + Msel*nu(i);
    }
    
    objoptII = ucminf_rcpp_splines(init_valII,y,se,opti,outofsample,parP10,nstates,d,k,W,restricted,hyper_tan,control);
    
    estpar = Rcpp::as<arma::vec>(objoptII["par"]);
    
    var_est_gamma = var(estpar(estpar.size()-k) + W.submat(0,1,W.n_rows-1,(k-1))*estpar.subvec(estpar.size()-k+1,estpar.size()-1));
    
    betatilde.row(h) = vecscal(estpar.subvec(0,11), var_est_gamma).t();
  }
  
  crit_min = vecscal(beta_hat, VARgamma) - colMeans(betatilde);
  opt_fun = as_scalar(sum(pow(crit_min,2)));
  
  return opt_fun;
}  




// [[Rcpp::export]]
Rcpp::List ucminf_rcpp_IndInf(const vec& init_val, const int& sim_h, const vec& beta_hat, const mat& eps_sim, double& VARgamma,
                              const bool& opti, const bool& outofsample, const double& parP10, int& d, const mat& W,
                              const bool& restricted, const bool& hyper_tan, const mat& se, const vec& states_noerr, const vec& init_val_CC,
                              const mat& Rsel, const vec& Msel, const vec& init_valII, Rcpp::List& control){
  
  // Extract R's optim function
  Rcpp::Environment stats("package:ucminf"); 
  Rcpp::Function ucminf = stats["ucminf"];
  
  // Call the optim function from R in C++ 
  Rcpp::List opt_results = ucminf(Rcpp::_["par"]    = init_val,
                                  // Make sure this function is not exported!
                                  Rcpp::_["fn"]     = Rcpp::InternalFunction(&IndInf_CC_rcpp),
                                  Rcpp::_["sim_h"] = sim_h,
                                  Rcpp::_["beta_hat"] = beta_hat,
                                  Rcpp::_["eps_sim"] = eps_sim,
                                  Rcpp::_["VARgamma"] = VARgamma,
                                  Rcpp::_["opti"] = opti,
                                  Rcpp::_["outofsample"] = outofsample,
                                  Rcpp::_["parP10"] = parP10,
                                  Rcpp::_["d"] = d,
                                  Rcpp::_["W"] = W,
                                  Rcpp::_["restricted"] = restricted,
                                  Rcpp::_["hyper_tan"] = hyper_tan,
                                  Rcpp::_["se"] = se,
                                  Rcpp::_["states_noerr"] = states_noerr,
                                  Rcpp::_["init_val_CC"] = init_val_CC,
                                  Rcpp::_["Rsel"] = Rsel,
                                  Rcpp::_["Msel"] = Msel,
                                  Rcpp::_["init_valII"] = init_valII,
                                  Rcpp::_["control"] = control);      // it does not take too many arguments!!
  
  // Extract out the estimated parameter values
  vec est_par = Rcpp::as<arma::vec>(opt_results[0]);      // estimated parameters
  double value = opt_results[1];      // value of the log-likelihood at the estimates
  
  // Return estimated values
  Rcpp::List ret;
  ret["par"] = est_par;
  ret["value"] = value;
  return ret;
}



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
  
  //Build T;
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
  
  //Bulild Z:
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
  
  //KF recursions
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

  //Bulild Z:
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

  //initialization
  for (int m = 0; m < draw_m; m++){
    gamma_t = as_scalar(exp(tau_hat(tau_hat.size()-1))*randn(1) + init_gamma);      // initialize gamma_t at the first observation of the cubic splines estimate

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


  // //loop over time
  for (int i = 0; i < len; i++){
    for (int m = 0; m < draw_m; m++){

      alpha_tm = vectorise(alpha_t.submat(0,m,alpha_t.n_rows-2,m));
      Pttm1_tm = list2mat(Pttm1_KF,m);
      gamma_tm = as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m));

      KF_results = KF_t_rcpp(tau_hat.subvec(0,tau_hat.size()-2),vectorise(y.col(i)),vectorise(se.row(i)),alpha_tm,Pttm1_tm,hyper_tan,gamma_tm,nstates);

      //compute weights
      p_den_m = KF_results["value"];
      p_density(m) = exp(-p_den_m);
      if (p_density(m) == 0){
        p_density(m) = 1e-20;
      }
      w_tilde(m) = w_t(m)*p_density(m);
    }

    //normalised weights
    w_t = w_tilde/sum(w_tilde);

    ESS(i) = as_scalar(1/(sum(pow(w_t,2))));
    CV(i) = sqrt(mean(pow(draw_m*w_t-1,2)));

    //resample
    w_t_power = pow(w_t,(log(sqrt(draw_m-1)/CV(i))));
    resample_set = stratifiedResampling_rcpp(w_t_power, draw_m);
    alpha_t = repcol(alpha_t,resample_set);
    Pttm1_KF = replist(Pttm1_KF,resample_set);

    w_t = ones(draw_m)/draw_m;

    //estimate state variable and variance
    for (int m = 0; m < draw_m; m++){
      att_BF(i) = att_BF(i) + w_t(m)*as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m));
      Ptt_BF(i) = Ptt_BF(i) + w_t(m)*pow(as_scalar(alpha_t.submat(alpha_t.n_rows-1,m,alpha_t.n_rows-1,m)),2);
    }
    Ptt_BF(i) = Ptt_BF(i) - pow(att_BF(i),2);

    if (i < len-1){
      //regenerate data
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

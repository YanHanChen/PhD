#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Dq0(double n, double f1, double f2, double g1,  double g2, double A) {
  double ans;
  if(f1+f2 == 0) {ans = 0;
  } else if ( ((g1*f2)/(2*f1) < g2) | (f1==0)){
    ans = ((n-1)/n)*(pow(g1,2)/(2*g2));
  } else{
    ans = ((n-1)/n)*(g1*(f1-1)/(2*(f2+1)));
  }
  return(ans);
}
// [[Rcpp::export]]
double Dq1_1(double n, double g1, double A) {
  double q1 = 0;
  double h2 = 0;
  if(A==1){
    h2 = 0;
  }else{
    for(int r = 1; r < n; r++){
      q1 = q1 + pow((1-A),r)/r;
    }
    h2 = (g1/n)*(pow(1-A,(-n+1)))*(-log(A)-q1);
  }
  return(h2);
}
// [[Rcpp::export]]
double Dq2(NumericMatrix tmpaL, double n, double t_bar) {
  double ans = 0;
  for(int i = 0; i < tmpaL.nrow(); i++){
    ans = ans + (tmpaL(i,2)*tmpaL(i,1)*tmpaL(i,0)*(tmpaL(i,0)-1)/n/(n-1));
  }
  //Rcout << "ans: " << ans;
  ans = pow(t_bar,2)/ans;
  return(ans);
}
// [[Rcpp::export]]
double delta(NumericMatrix del_tmpaL, double k, double n){
  double ans = 0;
  for(int i = 0;i < del_tmpaL.nrow(); i++){
    ans = ans +
      del_tmpaL(i,2)*del_tmpaL(i,1)*exp(Rf_lchoose(n-k-1,del_tmpaL(i,0)-1)-Rf_lchoose(n, del_tmpaL(i,0)));
  }
  return(ans);
}

// [[Rcpp::export]]
NumericVector RPD(NumericMatrix x , int n  , int m , NumericVector q) {
  int nrow = x.nrow();
  double tbar=0;
  NumericVector ghat(m);
  int qlength = q.length();
  NumericVector out(qlength);
  for (int i = 0; i < nrow; i++) {
    tbar += x(i, 0)*x(i, 1)/n;
  }
  //Rcpp::Rcout << "tbar in cpp: " << tbar << std::endl;
  for (int k = 0; k < m ; k++) {
    for (int i = 0; i < nrow; i++) {
      if ( x(i,0) >= k+1 && x(i,0) <= n-m+k+1 )
      {
        ghat[k] +=  x(i,1)*exp(Rf_lchoose(x(i,0), k+1)+Rf_lchoose(n-x(i,0), m-k-1)-Rf_lchoose(n, m)) ;
      }
      else
      {
        ghat[k] += 0 ;
      }
    }
  }
  //Rcpp::Rcout << "ghat: " << ghat << std::endl;
  for (int j = 0; j < qlength; j++ ){
    for(int k = 0; k < m; k++){
      if(q[j] == 0){
        out[j] = ghat[k] + out[j];
      }
      if(q[j] == 1){
        //Rcout << "q1 in cpp: " <<log ( (k+1) )<< std::endl;
        out[j] = -( (k+1) / (m*tbar) ) * log ( (k+1) / (m*tbar) ) * ghat[k] + out[j];
      }
      if(q[j] == 2){
        out[j] = pow( ( (k+1) / (m*tbar) ),2) * ghat[k] + out[j];
      }
    }
  }
  for(int j = 0; j < qlength; j++ ){
    if(q[j] == 0){
      out[j] = out[j] ;
    }
    if(q[j] == 1){
      out[j] = exp(out[j]);
    }
    if(q[j] == 2){
      out[j] = 1 / out[j];
    }
  }
  return out ;
}

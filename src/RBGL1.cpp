#include<RcppArmadillo.h>
#include<stdio.h>
#include<vector>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
//using namespace R;

// [[Rcpp::export()]]

Rcpp::List RBGL_1(arma::mat y, arma:: mat e, arma:: mat C, arma::mat g, arma:: mat w, arma:: vec z,int maxSteps, int n, int k,arma::vec hatBeta, arma:: mat hatEta, arma::vec hatAlpha, double hatTau, arma::vec hatV, arma::vec hatSg1,arma::vec hatSg2,arma::vec hatAta, arma::vec invSigAlpha0, double hatEtaSq1, double hatEtaSq2,double xi1, double xi2, double r1,double r2,double hatPhiSq,double a, double b, double alpha1,double gamma1, int progress)
{
  unsigned int q = e.n_cols,m = g.n_cols,p = w.n_cols, o = C.n_cols;
  arma::mat gsAlpha(maxSteps, q+o),
  gsBeta(maxSteps,m),
  gseta(maxSteps,p),
  gsAta(maxSteps,n),
  gsV(maxSteps, n*k),
  gsSg1(maxSteps, m),
  gsSg2(maxSteps, m)
    ;
  
  arma::vec gsEtaSq1(maxSteps),
  gsEtaSq2(maxSteps),
  gsTau(maxSteps),
  gsPhiSq(maxSteps)
    ;
  
  
  arma::mat mat0,mat1,mat2,mat3,ei,gi,wi;
  
  
  double meanAlpha;
  double varAlpha;
  
  double meanb;
  double varb;
  
  arma::vec muV, muS1;
  double muS2;
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2);
  arma::mat XgXgoV2(q,q);
  arma::rowvec RXgoV2(q);
  double RZoV;
  double tZZoV;
  
  for (int t = 0; t < maxSteps; t++) {
    
    
    mat0 = arma::repelem(e,k,1);
    mat1 = arma::repelem(g,k,1);
    mat2 = arma::repelem(w,k,1);
    
    
    // alpha|
    for(unsigned int j=0;j<q+o;j++){
      arma::vec res1, res11;
      double A0;
      A0 =0;
      double B0;
      B0=0;
      for(int i=0;i<n;i++){
        ei = mat0.rows((i*k),(i*k+k-1));
        ei.insert_cols(q, C);
        gi = mat1.rows((i*k),(i*k+k-1));
        wi = mat2.rows((i*k),(i*k+k-1));
        double tWWoV;
        tWWoV = arma::as_scalar((ei.col(j)/ hatV.subvec((i*k),(i*k+k-1))).t() * ei.col(j));
        A0 = A0+tWWoV;
        res1 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*arma::vectorise(hatEta)-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
        res11 = res1+ei.col(j)*hatAlpha(j);
        double RWoV;
        RWoV = arma::sum(ei.col(j) % (res11/ hatV.subvec((i*k),(i*k+k-1))));
        B0 = B0+ RWoV;
        
      }
      
      varAlpha = 1/(A0*hatTau/xi2Sq + invSigAlpha0(j));
      meanAlpha = varAlpha*B0*hatTau / xi2Sq;
      
      hatAlpha(j) = R::rnorm(meanAlpha,sqrt(varAlpha));
      
    }
    
    gsAlpha.row(t) = hatAlpha.t();
    
    // ata|
    arma::vec res;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(q, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));
      
      tZZoV = arma::as_scalar((z/hatV.subvec((i*k),(i*k+k-1))).t() * z);
      res = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*arma::vectorise(hatEta)-xi1*hatV.subvec((i*k),(i*k+k-1));
      RZoV = arma::sum(z% (res/hatV.subvec((i*k),(i*k+k-1))));
      double varAta;
      varAta= 1/(tZZoV*hatTau/xi2Sq+1/hatPhiSq);
      double meanAta;
      meanAta= varAta* RZoV * hatTau / xi2Sq;
      hatAta(i) = R::rnorm(meanAta, sqrt(varAta));
    }
    gsAta.row(t) = hatAta.t();
    
    //v|
    arma::vec resv;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(q, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));
      resv = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*arma::vectorise(hatEta)-z*hatAta(i);
      lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
      muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(resv));
      arma::vec v(k);
      for(int k0=0;k0<k;k0++){
        bool flag = true;
        while(flag){
          v(k0) = 1/rinvGauss(muV(k0), lambV);
          if(v(k0)<=0 || std::isinf(v(k0)) || std::isnan(v(k0))){
            if(progress != 0) Rcpp::Rcout << "v(k0) <= 0 or nan or inf" << std::endl;
            Rcpp::checkUserInterrupt();
          }else{
            flag = false;
          }
        }
      }
      hatV.subvec((i*k),(i*k+k-1)) = v;
    }
    
    
    gsV.row(t) = hatV.t();
    
    //s1|
    
    muS1 = std::sqrt(hatEtaSq1)/ arma::abs(hatBeta);
    for(unsigned int j = 0; j<m; j++){
      bool flag = true;
      while(flag){
        hatSg1(j) = 1/rinvGauss(muS1(j), hatEtaSq1);
        if(hatSg1(j)<=0 || std::isinf(hatSg1(j)) || std::isnan(hatSg1(j))){
          if(progress != 0) Rcpp::Rcout << "hatSg1(j): " << hatSg1(j) << std::endl;
          Rcpp::checkUserInterrupt();
        }else{
          flag = false;
        }
      }
    }
    gsSg1.row(t) = hatSg1.t();
    
    
    //s2|
    for(unsigned int j = 0; j<m; j++){
      muS2 = std::sqrt(hatEtaSq2)/ (arma::norm(hatEta.col(j)));
      bool flag = true;
      while(flag){
        hatSg2(j) = 1/rinvGauss(muS2, hatEtaSq2);
        if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
          if(progress != 0){
            Rcpp::Rcout << "hatSg2(j) = " << hatSg2(j) << " mu: " << muS2 << " lamb: " << hatEtaSq2 << std::endl;
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    gsSg2.row(t) = hatSg2.t();
    
    // Beta|
    
    for(unsigned int j=0;j<m;j++){
      arma::vec res2, res22;
      double A1;
      A1=0;
      double B1;
      B1=0;
      
      for(int i=0;i<n;i++){
        ei = mat0.rows((i*k),(i*k+k-1));
        ei.insert_cols(q, C);
        gi = mat1.rows((i*k),(i*k+k-1));
        wi = mat2.rows((i*k),(i*k+k-1));
        double XgXgoV1;
        XgXgoV1 = arma::as_scalar((gi.col(j)/ hatV.subvec((i*k),(i*k+k-1))).t() * gi.col(j));
        A1 = A1+XgXgoV1;
        res2 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*arma::vectorise(hatEta)-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
        res22 = res2+gi.col(j)*hatBeta(j);
        double RXgoV1;
        RXgoV1 = arma::sum(gi.col(j) % (res22/ hatV.subvec((i*k),(i*k+k-1))));
        B1 = B1+ RXgoV1;
        
      }
      
      varb = 1/(A1*hatTau/xi2Sq + 1/hatSg1(j));
      meanb = varb*B1*hatTau / xi2Sq;
      
      
      hatBeta(j)=R::rnorm(meanb,sqrt(varb));
      
    }
    
    gsBeta.row(t) = hatBeta.t();
    
    
    // eta|
    
    arma::vec meane;
    arma::mat varcove(q,q);
    
    arma:: mat Tau(q,q);
    Tau = Tau.eye();
    
    
    for(unsigned int j=0;j<m;j++){
      arma::mat A2(q,q);
      A2 = A2.zeros();
      arma:: vec B2;
      B2 = zeros<vec>(q);
      
      arma::vec res3, res33;
      
      for(int i=0;i<n;i++){
        ei = mat0.rows((i*k),(i*k+k-1));
        ei.insert_cols(q, C);
        gi = mat1.rows((i*k),(i*k+k-1));
        wi = mat2.rows((i*k),(i*k+k-1));
        XgXgoV2 = (wi.cols((j*q),(j*q+q-1)).each_col()/hatV.subvec((i*k),(i*k+k-1))).t()*wi.cols((j*q),(j*q+q-1));
        A2 = A2+XgXgoV2;
        res3 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*arma::vectorise(hatEta)-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
        res33 = res3+wi.cols((j*q),(j*q+q-1))*hatEta.col(j);
        RXgoV2 = arma::sum(wi.cols((j*q),(j*q+q-1)).each_col()%(res33/hatV.subvec((i*k),(i*k+k-1))),0);
        B2 = B2+RXgoV2.t();
        
      }
      
      
      varcove = arma::inv(1/hatSg2(j)*Tau+A2*hatTau/xi2Sq);
      meane = varcove*B2*hatTau/xi2Sq;
      
      
      hatEta.col(j) = mvrnormCpp(meane,varcove);
      
      
    }
    
    gseta.row(t) = arma::vectorise(hatEta).t();
    
    //etasq1;
    
    double shape2 = m+1;
    double rate2 = arma::accu(hatSg1)/2 + r1;
    hatEtaSq1 = R::rgamma(shape2, 1/rate2);
    gsEtaSq1(t) = hatEtaSq1;
    
    //etasq2;
    
    double shape21 = (m+m*q)/2+1;
    double rate21 = arma::accu(hatSg2)/2 + r2;
    hatEtaSq2 = R::rgamma(shape21, 1/rate21);
    gsEtaSq2(t) = hatEtaSq2;
    
    
    // phi.sq|
    double shapePhi, ratePhi;
    shapePhi = alpha1 + n/2;
    
    ratePhi = gamma1 + 0.5*(arma::accu(square(hatAta)));
    
    hatPhiSq = 1/R::rgamma(shapePhi, 1/ratePhi);
    
    gsPhiSq(t) = hatPhiSq;
    
    //tau|
    
    double shape = a + 3*n*k/2;
    double rest;
    arma::vec restt;
    rest = 0;
    double f;
    f=0;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(q, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));
      restt = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*arma::vectorise(hatEta)-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
      double ResSqoV;
      ResSqoV = arma::accu(arma::square(restt)/hatV.subvec((i*k),(i*k+k-1)));
      rest = rest+ResSqoV/(2*xi2Sq);
      f = f+ arma::accu(hatV.subvec((i*k),(i*k+k-1)));
    }
    double rate = b + f + rest/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(t) = hatTau;
    
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("GS.alpha") = gsAlpha,
    Rcpp::Named("GS.beta") = gsBeta,
    Rcpp::Named("GS.eta") = gseta,
    Rcpp::Named("GS.ata") = gsAta,
    Rcpp::Named("GS.v") = gsV,
    Rcpp::Named("GS.s1") = gsSg1,
    Rcpp::Named("GS.s2") = gsSg2,
    Rcpp::Named("GS.eta21.sq") = gsEtaSq1,
    Rcpp::Named("GS.eta22.sq") = gsEtaSq2,
    Rcpp::Named("GS.phi.sq") = gsPhiSq,
    Rcpp::Named("GS.tau") = gsTau
  
  );
  
}
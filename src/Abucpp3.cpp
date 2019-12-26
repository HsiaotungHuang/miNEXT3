#include <Rcpp.h>


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

// [[Rcpp::export]]
double Hypergeometric(int K, int k, int N,  int n) {
  return exp(Rf_lchoose(K,k)+Rf_lchoose(N-K,n-k)-Rf_lchoose(N,n));
  //return Rf_choose(K,k)*Rf_choose(N-K,n-k)/Rf_choose(N,n);
}

// [[Rcpp::export]]
double fk_3(int k1, int k2,int k3, int m1, int m2,int m3, NumericVector x1, NumericVector y1, NumericVector z1){
  int n1 = sum(x1);
  int n2 = sum(y1);
  int n3 = sum(z1);
  NumericVector x = x1[(x1>=k1) & (y1>=k2)& (z1>=k3)];
  NumericVector y = y1[(x1>=k1) & (y1>=k2)& (z1>=k3)];
  NumericVector z = z1[(x1>=k1) & (y1>=k2)& (z1>=k3)];
  
  double output = 0;
  if(sum((x1>=k1) & (y1>=k2)& (z1>=k3))==0){
    output = output;
  }else{
    for(int i = 0; i<x.size(); i++){
      output = output + Hypergeometric(x[i], k1, n1, m1)*Hypergeometric(y[i], k2, n2, m2)*Hypergeometric(z[i], k3, n3, m3);
    }
  }
  return output;
}

// [[Rcpp::export]]
NumericVector D_share_yhc(NumericVector xi,NumericVector yi,NumericVector zi,double m1, double m2, double m3,NumericVector q){
  NumericVector fk123(q.size());
  double max1 = max(xi);
  double max2 = max(yi);
  double max3 = max(zi);
  double m1loop = std::min(m1,max1);
  double m2loop = std::min(m2,max2);
  double m3loop = std::min(m3,max3);
  double tmp = 0;
  NumericVector output(q.size());
  
  for(int k1 = 0; k1<(m1loop+1) ; k1++){
    for(int k2 = 0; k2<(m2loop+1) ; k2++){
      for(int k3 = 0; k3<(m3loop+1) ; k3++){
        if((k1==0) & ((k2 == 0) & (k3 == 0))){
          fk123 = fk123;
        }else{
          tmp = fk_3(k1, k2, k3, m1, m2, m3, xi, yi, zi);
          for(int qi = 0; qi< q.size(); qi++){
            if(q[qi]==1){ 
              fk123[qi] = fk123[qi] - ((k1+k2+k3)/(m1+m2+m3))*log((k1+k2+k3)/(m1+m2+m3))*tmp;
            }else{
              fk123[qi] = fk123[qi] + pow(((k1+k2+k3)/(m1+m2+m3)), q[qi])*tmp;
            }
          }
        }
      }
    }
  }
  for(int qi = 0; qi< q.size(); qi++){
    if(q[qi]==1){ 
      output[qi] = exp(fk123[qi]);
    }else{
      output[qi] = pow(fk123[qi], 1/(1-q[qi]));
    }
  }
  return output;
}
// [[Rcpp::export]]
double D0_rare_yhc(NumericVector xi,NumericVector yi,NumericVector zi,double m1, double m2, double m3){
  double output = 0;
  double n1 = sum(xi);
  double n2 = sum(yi);
  double n3 = sum(zi);
  
  for(int i = 0; i < xi.size(); i++){
    output = output + 1 - exp(Rf_lchoose(n1 - xi[i], m1)- Rf_lchoose(n1,m1))*
      exp(Rf_lchoose(n2 - yi[i], m2)- Rf_lchoose(n2,m2))*exp(Rf_lchoose(n3 - zi[i], m3)- Rf_lchoose(n3,m3));
  }
  return output;
}

// [[Rcpp::export]]
double D_q0_in_3(NumericVector xi,NumericVector yi,NumericVector zi,double m1, double m2, double m3){
  double output = 0;
  double n1 = sum(xi);
  double n2 = sum(yi);
  double n3 = sum(zi);
  
  for(int i = 0; i < xi.size(); i++){
    output = output + 1 - exp(Rf_lchoose(n1 - xi[i], m1)- Rf_lchoose(n1,m1))*
      exp(Rf_lchoose(n2 - yi[i], m2)- Rf_lchoose(n2,m2))*exp(Rf_lchoose(n3 - zi[i], m3)- Rf_lchoose(n3,m3));
  }
  return output;
}


// [[Rcpp::export]]
double q2_p_cpp(NumericVector x1, NumericVector y1, NumericVector z1,double m1, double m2, double m3,double n1, double n2, double n3){
  double p1=0;double p2=0;double p3=0;
  double p12=0;double p23=0;double p13=0;
  NumericVector outp(6);
  
  for(int i = 0; i<x1.size() ; i++){
    p1=p1+x1[i]*(x1[i]-1)/(n1*(n1-1));
  }
  for(int i = 0; i<y1.size() ; i++){
    p2=p2+y1[i]*(y1[i]-1)/(n2*(n2-1));
  }
  for(int i = 0; i<z1.size() ; i++){
    p3=p3+z1[i]*(z1[i]-1)/(n3*(n3-1));
  }
  
  for(int i = 0; i<x1.size() ; i++){
    p12=p12+x1[i]*y1[i]/(n1*n2);
  }
  
  for(int i = 0; i<x1.size() ; i++){
    p13=p13+x1[i]*z1[i]/(n1*n3);
  }
  for(int i = 0; i<y1.size() ; i++){
    p23=p23+y1[i]*z1[i]/(n2*n3);
  }
  outp[0]=p1;outp[1]=p2;outp[2]=p3;
  outp[3]=p12;outp[4]=p23;outp[5]=p13;
  
  return outp(6);
}

// [[Rcpp::export]]
double D_q01_in_3(NumericVector xi,NumericVector yi,NumericVector zi,double m1, double m2, double m3,double q){
  //NumericVector xi = X(_,0);
  //NumericVector yi = X(_,1);
  double fk123 = 0;
  double output = 0;
  
  double max1 = max(xi);
  double max2 = max(yi);
  double max3 = max(zi);
  double m1loop = std::min(m1,max1);
  double m2loop = std::min(m2,max2);
  double m3loop = std::min(m3,max3);
  
  
  if(q==1){
    for(int k1 = 0; k1<(m1loop+1) ; k1++){
      for(int k2 = 0; k2<(m2loop+1) ; k2++){
        for(int k3 = 0; k3<(m3loop+1) ; k3++){
          if((k1 == 0) & ((k2 == 0) & (k3 == 0))){
            fk123 = fk123; 
          }else{
            fk123 = fk123 -((k1+k2+k3)/(m1+m2+m3))*log((k1+k2+k3)/(m1+m2+m3))*fk_3(k1, k2, k3, m1, m2, m3, xi, yi, zi);
          }
        }
      }
    }
    output = exp(fk123);
    if ((m1loop==0) & ((m2loop==0) & (m3loop==0))) output=1;
  }else{
    for(int k1 = 0; k1<(m1loop+1) ; k1++){
      for(int k2 = 0; k2<(m2loop+1) ; k2++){
        for(int k3 = 0; k3<(m3loop+1); k3++){
          if((k1 == 0) & ((k2 == 0) & (k3 == 0))){
            fk123 = fk123; 
          }else{
            fk123 = fk123 + pow(((k1+k2+k3)/(m1+m2+m3)), q)*fk_3(k1, k2, k3, m1, m2, m3, xi, yi, zi);
          }
        }
      }
    }
    output = pow(fk123, 1/(1-q));
    if ((m1loop==0) & ((m2loop==0) & (m3loop==0))) output=0;
  }
  return output;
}


double h0_3_1cpp(double pi1, double pi2, double pi3,int m1,int m2,int m3s, int n3 ){
  double output = 0;
  output = pow((1-pi1),m1)*pow((1-pi2),m2)*pow((1-pi3),n3)*(1-pow((1-pi3),m3s));
  return output;
}


double h0_3_2cpp(double pi1, double pi2, double pi3,int m1,int m2s,int n2 ,int m3s, int n3 ){
  double output = 0;
  output = pow((1-pi1),m1)*pow((1-pi2),n2)*pow((1-pi3),n3)*(1-pow((1-pi2),m2s)*pow((1-pi3),m3s));
  return output;
}


// [[Rcpp::export]]
double h0_3_1hat_cpp(NumericVector pi1, NumericVector pi2, NumericVector pi3, int m1, int m2,int m3s, int n1, int n2, int n3){
  double output_all= 0; 
  //  double output_sh = 0;
  if(m1>=0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h0_3_1cpp(pi1_tmp[i],pi2_tmp[i],pi3_tmp[i],m1,m2,m3s,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3>0)];
    double sum011 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum011 = sum011 +  h0_3_1cpp(0,pi2_tmp[i],pi3_tmp[i],m1,m2,m3s,n3)/((1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3==0)];
    double sum110 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum110 = sum110 +  h0_3_1cpp(pi1_tmp[i],pi2_tmp[i],0,m1,m2,m3s,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3>0)];
    double sum101 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum101 = sum101 +  h0_3_1cpp(pi1_tmp[i],0,pi3_tmp[i],m1,m2,m3s,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2==0)& (pi3>0)];
    double sum001 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum001 = sum001 +  h0_3_1cpp(0,0,pi3_tmp[i],m1,m2,m3s,n3)/(1-pow(1-pi3_tmp[i], n3));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3==0)];
    double sum010 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum010 = sum010 +  h0_3_1cpp(0,pi3_tmp[i],0,m1,m2,m3s,n3)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3==0)];
    double sum100 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum100 = sum100 +  h0_3_1cpp(pi1_tmp[i],0,0,m1,m2,m3s,n3)/(1-pow(1-pi1_tmp[i], n1));
    }
    
    output_all =sumsh+sum011+sum110+sum101+sum001+sum010+sum100;
    
    //    output_sh = sumsh;
  }
 // else if(m1 == 0 & m2==0){
    //###call iNEXT
 // }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}


// [[Rcpp::export]]
double h0_3_2hat_cpp(NumericVector pi1, NumericVector pi2, NumericVector pi3, int m1, int m2s,int m3s, int n1, int n2, int n3){
  double output_all= 0; 
  //  double output_sh = 0;
  if(m1>=0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h0_3_2cpp(pi1_tmp[i],pi2_tmp[i],pi3_tmp[i],m1,m2s,n2,m3s,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3>0)];
    double sum011 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum011 = sum011 +  h0_3_2cpp(0,pi2_tmp[i],pi3_tmp[i],m1,m2s,n2,m3s,n3)/((1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3==0)];
    double sum110 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum110 = sum110 +  h0_3_2cpp(pi1_tmp[i],pi2_tmp[i],0,m1,m2s,n2,m3s,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3>0)];
    double sum101 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum101 = sum101 +  h0_3_2cpp(pi1_tmp[i],0,pi3_tmp[i],m1,m2s,n2,m3s,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2==0)& (pi3>0)];
    double sum001 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum001 = sum001 +  h0_3_2cpp(0,0,pi3_tmp[i],m1,m2s,n2,m3s,n3)/(1-pow(1-pi3_tmp[i], n3));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3==0)];
    double sum010 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum010 = sum010 +  h0_3_2cpp(0,pi3_tmp[i],0,m1,m2s,n2,m3s,n3)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3==0)];
    double sum100 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum101 = sum100 +  h0_3_2cpp(pi1_tmp[i],0,0,m1,m2s,n2,m3s,n3)/(1-pow(1-pi1_tmp[i], n1));
    }
    
    output_all =sumsh+sum011+sum010+sum100+sum001+sum010+sum100;
    
    //    output_sh = sumsh;
  }
 // else if(m1 == 0 & m2s==0){
    //###call iNEXT
 // }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}

// [[Rcpp::export]]
double Efk_q1_3(double pi1, double pi2, double pi3,int m1,int m2,int m3,int k1,int k2,int k3){
  double result;
  double r1;double r2;double r3;
  if((k1 == 0) & (k2==0) & (k3==0)){
    result = 0;
  }
  else{
    if((pi1 == 0) & (k1==0)){r1=1;}
    else{
      r1=Rf_dbinom(k1, m1 , pi1, 0);
    }
    
    if((pi2 == 0) & (k2==0)){r2=1;}
    else{
      r2=Rf_dbinom(k2, m2 , pi2, 0);
    }
    
    if((pi3 == 0) & (k3==0)){r3=1;}
    else{
      r3=Rf_dbinom(k3, m3 , pi3, 0);
    }
    result = r1*r2*r3;
  }
  return result;
}

// [[Rcpp::export]]
double h1_3_1cpp(double pi1, double pi2, double pi3,double m1,double m2,double m3,double n1,double n2,double n3){
  
  double tmp1 = 0;
  double tmp2 = 0;
  for(int k3=0; k3 <= m3; k3++){
    for(int k2=0; k2 <= m2; k2++){
      for(int k1=0; k1 <= m1; k1++){
        if((k1 == 0) & (k2 == 0)& (k3 == 0)){
          tmp1 = 0;
          }else{ 
          tmp1 = tmp1 + (k1+k2+k3)/(m1+m2+m3)*log((k1+k2+k3)/(m1+m2+m3))*Efk_q1_3(pi1,pi2,pi3,m1,m2,m3,k1,k2,k3); }
      }
    }
  }
  
  tmp1 = -tmp1;
  for(int k3=0; k3 <= n3; k3++){
    for(int k2=0; k2 <= m2; k2++){
      for(int k1=0; k1 <= m1; k1++){
        if((k1 == 0) & (k2 == 0)& (k3 == 0)){
          tmp2 = 0;
          }else{ 
          tmp2 = tmp2 + (k1+k2+k3)/(m1+m2+m3)*log((k1+k2+k3)/(m1+m2+m3))*Efk_q1_3(pi1,pi2,pi3,m1,m2,n3,k1,k2,k3); }
      }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}

// [[Rcpp::export]]
double h1_3_2cpp(double pi1, double pi2, double pi3,double m1,double m2,double m3,double n1,double n2,double n3){
  
  double tmp1 = 0;
  double tmp2 = 0;
  for(int k3=0; k3 <= m3; k3++){
    for(int k2=0; k2 <= m2; k2++){
      for(int k1=0; k1 <= m1; k1++){
        if((k1 == 0) & (k2 == 0)& (k3 == 0)){ tmp1 = 0; }
        else{ tmp1 = tmp1 + (k1+k2+k3)/(m1+m2+m3)*log((k1+k2+k3)/(m1+m2+m3))*Efk_q1_3(pi1,pi2,pi3,m1,m2,m3,k1,k2,k3); }
      }
    }
  }
  
  tmp1 = -tmp1;
  
  for(int k3=0; k3 <= n3; k3++){
    for(int k2=0; k2 <= n2; k2++){
      for(int k1=0; k1 <= m1; k1++){
        if((k1 == 0) & (k2 == 0)& (k3 == 0)){ tmp2 = 0; }
        else{ tmp2 = tmp2 + (k1+k2+k3)/(m1+n2+n3)*log((k1+k2+k3)/(m1+n2+n3))*Efk_q1_3(pi1,pi2,pi3,m1,n2,n3,k1,k2,k3); }
      }
    }
  }
  tmp2 = -tmp2;
  // Rcout << "h1 function, tmp1= " << tmp1 << std::endl;
  // Rcout << "h1 function, tmp2= " << tmp2 << std::endl;
  double result = tmp1-tmp2;
  return result;
}


// [[Rcpp::export]]
double h1_3_1hat_cpp(NumericVector pi1, NumericVector pi2, NumericVector pi3, double m1, double m2,double m3, double n1, double n2, double n3){
  double output_all= 0; 
  //int m3s=m3-n3;
  //  double output_sh = 0;
  if(m1>=0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1_3_1cpp(pi1_tmp[i],pi2_tmp[i],pi3_tmp[i],m1,m2,m3,n1,n2,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3>0)];
    double sum011 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum011 = sum011 +  h1_3_1cpp(0,pi2_tmp[i],pi3_tmp[i],m1,m2,m3,n1,n2,n3)/((1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3==0)];
    double sum110 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum110 = sum110 +  h1_3_1cpp(pi1_tmp[i],pi2_tmp[i],0,m1,m2,m3,n1,n2,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3>0)];
    double sum101 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum101 = sum101 +  h1_3_1cpp(pi1_tmp[i],0,pi3_tmp[i],m1,m2,m3,n1,n2,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2==0)& (pi3>0)];
    double sum001 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum001 = sum001 +  h1_3_1cpp(0,0,pi3_tmp[i],m1,m2,m3,n1,n2,n3)/(1-pow(1-pi3_tmp[i], n3));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3==0)];
    double sum010 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum010 = sum010 +  h1_3_1cpp(0,pi3_tmp[i],0,m1,m2,m3,n1,n2,n3)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3==0)];
    double sum100 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum100 = sum100 +  h1_3_1cpp(pi1_tmp[i],0,0,m1,m2,m3,n1,n2,n3)/(1-pow(1-pi1_tmp[i], n1));
    }
    output_all =sumsh+sum011+sum010+sum100+sum001+sum010+sum100;
    
    //    output_sh = sumsh;
  }
  // else if(m1 == 0 & m2 == 0 ){
  //   //###call iNEXT
  // }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}




// [[Rcpp::export]]
double h1_3_2hat_cpp(NumericVector pi1, NumericVector pi2, NumericVector pi3, double m1, double m2,double m3, double n1, double n2, double n3){
  double output_all= 0; 
  //int m3s=m3-n3;
  //  double output_sh = 0;
  if(m1>=0){
    NumericVector pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3>0)];
    NumericVector pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3>0)];
    double sumsh = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sumsh = sumsh +  h1_3_2cpp(pi1_tmp[i],pi2_tmp[i],pi3_tmp[i],m1,m2,m3,n1,n2,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3>0)];
    double sum011 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum011 = sum011 +  h1_3_2cpp(0,pi2_tmp[i],pi3_tmp[i],m1,m2,m3,n1,n2,n3)/((1-pow(1-pi2_tmp[i], n2))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2>0)& (pi3==0)];
    double sum110 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum110 = sum110 +  h1_3_2cpp(pi1_tmp[i],pi2_tmp[i],0,m1,m2,m3,n1,n2,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi2_tmp[i], n2)));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3>0)];
    double sum101 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum101 = sum101 +  h1_3_2cpp(pi1_tmp[i],0,pi3_tmp[i],m1,m2,m3,n1,n2,n3)/((1-pow(1-pi1_tmp[i], n1))*(1-pow(1-pi3_tmp[i], n3)));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2==0)& (pi3>0)];
    pi2_tmp = pi2[(pi1==0) & (pi2==0)& (pi3>0)];
    pi3_tmp = pi3[(pi1==0) & (pi2==0)& (pi3>0)];
    double sum001 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum001 = sum001 +  h1_3_2cpp(0,0,pi3_tmp[i],m1,m2,m3,n1,n2,n3)/(1-pow(1-pi3_tmp[i], n3));
    }
    
    pi1_tmp = pi1[(pi1==0) & (pi2>0)& (pi3==0)];
    pi2_tmp = pi2[(pi1==0) & (pi2>0)& (pi3==0)];
    pi3_tmp = pi3[(pi1==0) & (pi2>0)& (pi3==0)];
    double sum010 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum010 = sum010 +  h1_3_2cpp(0,pi3_tmp[i],0,m1,m2,m3,n1,n2,n3)/(1-pow(1-pi2_tmp[i], n2));
    }
    
    pi1_tmp = pi1[(pi1>0) & (pi2==0)& (pi3==0)];
    pi2_tmp = pi2[(pi1>0) & (pi2==0)& (pi3==0)];
    pi3_tmp = pi3[(pi1>0) & (pi2==0)& (pi3==0)];
    double sum100 = 0;
    for(int i=0; i < pi1_tmp.size(); i++){
      sum100 = sum100 +  h1_3_2cpp(pi1_tmp[i],0,0,m1,m2,m3,n1,n2,n3)/(1-pow(1-pi1_tmp[i], n1));
    }
    
    output_all =sumsh+sum011+sum010+sum100+sum001+sum010+sum100;
    
    //    output_sh = sumsh;
  }
  // else if(m1 == 0 & m2 == 0 ){
  //   //###call iNEXT
  // }
  //  NumericVector output = NumericVector::create(output_all, output_sh);
  double output = output_all;
  return output;
}



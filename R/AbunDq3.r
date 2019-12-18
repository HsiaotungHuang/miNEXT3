


#######################################################################################################################
#' R code for q=0,1, 3-habitat mixture of rarefactioncurves
#' @param x1 a vector/matrix/list of species abundance frequency
#' @param mm an int matrix for number of individuals for x1,x2,x3 
#' @param q  a numerical value for Dq
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution
#' @examples
#' data(Abudata)
#' data1<-Abudata$data
#' mm1<-Abudata$mm
#' Dq_in(data1[,1],data1[,2],data1[,3],mm1,0)
#' @export
Dq_in<-function(x1,x2,x3,mm,q){
  out = apply(mm,1,function(mm) {
    mm1 = mm[1]
    mm2 = mm[2]
    mm3 = mm[3]
    
    print(paste("mm1=",mm1,"mm2=",mm2,"mm3=",mm3,"q=",q))
    D_share(x1, x2,x3, mm1, mm2,mm3, q)
  })
  out
}


#######################################################################################################################
#' R code for q=0,1  3-habitat mixture of rarefaction curves
#' @param x1 a vector/matrix/list of species abundance frequency
#' @param mm an int matrix for number of individuals for x1,x2,x3 
#' @param q  a numerical value for Dq
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution
#' @examples
#' data(Abudata)
#' data1<-Abudata$data
#' mm1<-Abudata$mm
#' Dq_in_yhc(data1[,1],data1[,2],data1[,3],mm1,0)
#' @export
Dq_in_yhc<-function(x1,x2,x3,mm,q){
  q1pls <- q[!c(q==0|q==2)]
  out1pls <- apply(mm,1,function(i) {
    D_share_yhc(x1, x2,x3, i[1], i[2],i[3],q1pls)
  })
  if(length(q1pls)>=2) out1pls <- t(out1pls)
  out0 <- apply(mm,1,function(i) {
    D0_rare_yhc(x1, x2,x3, i[1], i[2],i[3])
  })
  out2 <- D2_yhc(data = cbind(x1,x2,x3),mm)
  cbind(out0,out1pls,out2)
}

#######################################################################################################################
#' R code for q=2, 3-habitat mixture of rarefaction/extrapolation curves
#' @import dplyr
#' @param x1 a vector/matrix/list of species abundance frequency
#' @param mm an int matrix for number of individuals for x1,x2,x3 
#' @param q  a numerical value for Dq
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution
#' @examples
#' data(Abudata)
#' data1<-Abudata$data
#' mm1<-Abudata$mm
#' Dq_in_yhc(data1[,1],data1[,2],data1[,3],mm1,0)
#' @export
D2_yhc <- function(data,mm){
  N <- ncol(data)
  n <- colSums(data)
  ch2 <- combn(x = 1:N,m = 2,simplify = T) %>% t()
  pk <- sapply(1:N, function(j){
    sum(data[,j]/n[j]*(data[,j]-1)/(n[j]-1))
  })
  pkj <- apply(ch2, 1, function(j){
    sum(data[,j[1]]/n[j[1]] * data[,j[2]]/n[j[2]])
  })
  
  subfun <- function(mm){
    tmp1 <- sum(mm*(mm-1)*pk/(sum(mm)^2))
    tmp2 <- 2*sapply(1:nrow(ch2),function(k){
      mm[ch2[k,]] %>% prod(.) * pkj[k]
    }) %>% sum(.)/(sum(mm)^2)
    return((1/sum(mm) + tmp1 + tmp2)^(-1))
  }
  apply(mm, 1, subfun)
}

D_q0_ext1<-function(x1,x2,x3,p1,p2,p3,mm,n,q){
  D_q0_ext1_value<-0
  n1<-n[1]
  n2<-n[2]
  n3<-n[3]
  mm3<-mm[,3]
  mm3<-mm3[mm3>n3]
  if (length(mm3)==0 | length(mm)<length(mm3)) {stop("function D0_ext1 :all m3 should be larger than n3")}
  else{
    mm.1<-cbind(mm[,1],mm[,2],n3)
    q0_tmp1<-apply(mm.1,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3 =mm[3]
      out<-D_q01_in_3(p1,p2,p3, mm1, mm2, mm3,0)
    })
    mmext<-cbind(mm[,1],mm[,2],mm[,3]-n3)
    h0_3_1hat = apply(mmext,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3s =mm[3]
      out<-h0_3_1hat_cpp(p1,p2,p3, mm1, mm2, mm3s,n1,n2,n3)
    })
    
    D_q0_ext1_value<-q0_tmp1+h0_3_1hat
  }
    return(D_q0_ext1_value)
}

##input those mm that m1<=n1, but m2>n2, m3>n3
D_q0_ext2<-function(x1,x2,x3,p1,p2,p3,mm,n,q){
  D_q0_ext2_value<-0
  n1<-n[1]
  n2<-n[2]
  n3<-n[3]
  mm3<-mm[,3]
  mm3<-mm3[mm3>n3]
  mm2<-mm[,2]
  mm2<-mm2[mm2>n3]
  if (length(mm3)==0 | length(mm)<length(mm3)) stop("function D0_ext2 :all m3 should be larger than n3")
  else if (length(mm2)==0 | length(mm)<length(mm2)) stop("function D0_ext2 :all m2 should be larger than n2")
  else{
    mm.1<-cbind(mm[,1],n2,n3)
    
    q0_tmp1<-apply(mm.1,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3 =mm[3]
      out<-D_q01_in_3(p1,p2,p3, mm1, mm2, mm3,0)
    })
    
    mmext<-cbind(mm[,1],mm[,2]-n2,mm[,3]-n3)
    h0_3_2hat = apply(mmext,1,function(mm) {
      mm1 = mm[1]
      mm2s = mm[2]
      mm3s =mm[3]
      out<-h0_3_2hat_cpp(p1,p2,p3, mm1, mm2s, mm3s,n1,n2,n3)
    })
    
    D_q0_ext2_value<-q0_tmp1+h0_3_2hat
  }
  return(D_q0_ext2_value)
}


D_q1_ext1<-function(x1,x2,x3,p1,p2,p3,mm,n,q){
  D_q1_ext1_value<-0
  n1<-n[1]
  n2<-n[2]
  n3<-n[3]
  mm3<-mm[,3]
  mm3<-mm3[mm3>n3]
  if (length(mm3)==0 | length(mm)<length(mm3)) stop("function D0_ext1 :all m3 should be larger than n3")
  else{
    mm.1<-cbind(mm[,1],mm[,2],n3)
    q1_tmp1<-apply(mm.1,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3 =mm[3]
      out<-D_q01_in_3(p1,p2,p3, mm1, mm2, mm3,1)
    })
    
    
    h1_3_1hat = apply(mmext,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3 =mm[3]
      out<-h1_3_1hat_cpp(p1,p2,p3, mm1, mm2, mm3,n1,n2,n3)
    })
    
    D_q1_ext1_value<-q1_tmp1+exp(h1_3_1hat)
  }
  return(D_q1_ext1_value)
}


D_q1_ext2<-function(x1,x2,x3,p1,p2,p3,mm,n,q){
  D_q1_ext2_value<-0
  n1<-n[1]
  n2<-n[2]
  n3<-n[3]
  mm3<-mm[,3]
  mm3<-mm3[mm3>n3]
  if (length(mm3)==0 | length(mm)<length(mm3)) stop("function D0_ext1 :all m3 should be larger than n3")
  else{
    mm.1<-cbind(mm[,1],mm[,2],n3)
    
    q1_tmp1<-apply(mm.1,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3 =mm[3]
      out<-D_q01_in_3(p1,p2,p3, mm1, mm2, mm3,1)
    })
    
    h1_3_2hat = apply(mmext,1,function(mm) {
      mm1 = mm[1]
      mm2 = mm[2]
      mm3 =mm[3]
      out<-h1_3_2hat_cpp(p1,p2,p3, mm1, mm2, mm3,n1,n2,n3)
    })
    
    D_q1_ext2_value<-q1_tmp1+exp(h1_3_2hat)
  }
  return(D_q1_ext2_value)
}





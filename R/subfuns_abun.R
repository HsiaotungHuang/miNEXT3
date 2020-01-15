#' @import dplyr
#' @importFrom chaoUtility Boot_p
mix3 <- function(data,ms){
  n <- colSums(data)
  if( sum( apply(ms, 1, function(m) sum(m<=n)) ==3) < nrow(ms)){
    datap <- Boot_p(x = data,Bootype = "JADE",datatype = "abundance") %>% sapply(., function(x){
      x[1:nrow(data)]
    })
  }else{
    datap <- NULL
  }
  apply(ms,1,function(m){
    mix3_each(data,m,n,datap)
  }) %>% t()
}

mix3_each <- function(data,m,n,datap=NULL){
  if(sum(m>n)==0){
    q0 <- D0_rare(data[,1], data[,2],data[,3], m[1], m[2],m[3])
    q1 <- D_rare(data[,1], data[,2],data[,3], m[1], m[2],m[3],1)
    q2 <- D2(data,m)
  }else if(sum(m[-1]>n[-1])==1){
    if(m[2] > n[2]){
      data <- data[,c(1,3,2)]
      datap <- datap[,c(1,3,2)]
      n <- n[c(1,3,2)]
      m <- m[c(1,3,2)]
    }
    q0 <- D_q0_ext1(data,datap,m,n)
    q1 <- D_q1_ext1(data,datap,m,n)
    q2 <- D2(data,m)
  }else{
    q0 <- D_q0_ext2(data,datap,m,n)
    q1 <- D_q1_ext2(data,datap,m,n)
    q2 <- D2(data,m)
  }
  return(c(q0,q1,q2))
}

# q = 2, work for both interpolation and extrapolation
#' @importFrom utils combn
#' @import dplyr
D2 <- function(data,mm){
  N <- ncol(data)
  n <- colSums(data)
  ch2 <- combn(x = 1:N,m = 2,simplify = T) %>% t()
  pk <- sapply(1:N, function(j){
    sum(data[,j]/n[j]*(data[,j]-1)/(n[j]-1))
  })
  pkj <- apply(ch2, 1, function(j){
    sum(data[,j[1]]/n[j[1]] * data[,j[2]]/n[j[2]])
  })
  
  tmp1 <- sum(mm*(mm-1)*pk/(sum(mm)^2))
  tmp2 <- 2*sapply(1:nrow(ch2),function(k){
    mm[ch2[k,]] %>% prod(.) * pkj[k]}
  ) %>% sum(.)/(sum(mm)^2)
  return((1/sum(mm) + tmp1 + tmp2)^(-1))
  
  # apply(mm, 1, subfun)
}

# q = 0, for two rarefaction and one extrapolation
D_q0_ext1<-function(data,datap,m,n){
  q0_tmp1 <- D0_rare(data[,1],data[,2],data[,3], m[1], m[2], n[3])
  h0_3_1hat <- h0_3_1hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2], m[3]-n[3],n[1],n[2],n[3])
  
  D_q0_ext1_value<-q0_tmp1+h0_3_1hat
  
  return(D_q0_ext1_value) 
}
# q = 0, for one rarefaction and two extrapolation
D_q0_ext2<-function(data,datap,m,n){
  q0_tmp1<-D0_rare(data[,1],data[,2],data[,3], m[1], n[2], n[3])
  h0_3_2hat <- h0_3_2hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2]-n[2], m[3]-n[3],n[1],n[2],n[3])
  D_q0_ext2_value<-q0_tmp1+h0_3_2hat
  return(D_q0_ext2_value)
}
# q = 1, for two rarefaction and one extrapolation
D_q1_ext1<-function(data,datap,m,n){
  q1_tmp1 <- D_rare(data[,1],data[,2],data[,3], m[1], m[2], n[3],1)
  h1_3_1hat <- h1_3_1hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2], m[3],n[1],n[2],n[3])
  
  D_q1_ext1_value<-q1_tmp1*exp(h1_3_1hat)
  return(D_q1_ext1_value)
}
# q = 1, for one rarefaction and two extrapolation
D_q1_ext2<-function(data,datap,m,n){
  q1_tmp1 <- D_rare(data[,1],data[,2],data[,3], m[1], n[2], n[3],1)
  h1_3_2hat <- h1_3_2hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2], m[3],n[1],n[2],n[3])
  D_q1_ext2_value<-q1_tmp1*exp(h1_3_2hat)
  return(D_q1_ext2_value)
}

# Functions for bootstrap
#' @importFrom chaoUtility Boot_p
#' @importFrom stats rmultinom
Abun_CreatBootstrapSample <- function(data, nboots = 0){
  data <- data[rowSums(data)>0,]
  if(nboots>1){
    n <- colSums(data)
    x1 = data[, 1]
    x2 = data[, 2]
    x3 = data[, 3]
    
    p1_est = Boot_p(x = x1,datatype = "abundance")
    p2_est = Boot_p(x = x2,datatype = "abundance")
    p3_est = Boot_p(x = x3,datatype = "abundance")
    data_p = data
    data_p[,1]<-p1_est[1:nrow(data_p)]
    data_p[,2]<-p2_est[1:nrow(data_p)]
    data_p[,3]<-p3_est[1:nrow(data_p)]
    undetec1 <- p1_est[-c(1:nrow(data_p))]
    undetec2 <- p2_est[-c(1:nrow(data_p))]
    undetec3 <- p3_est[-c(1:nrow(data_p))]
    
    
    f1 = sum(rowSums(data)==1)
    f2 = sum(rowSums(data)==2)
    f0hat<-Chat1_f0Fun(f1,f2,sum(n))[2]
    
    data_p = matrix(data = 0,nrow = f0hat, ncol = 3,
                    dimnames = list(NULL, names(data_p))) %>% rbind(data_p,.)
    zero1 <- which(data_p[,1]==0)
    zero2 <- which(data_p[,2]==0)
    zero3 <- which(data_p[,3]==0)
    
    data_boot <- lapply(1:nboots,function(B){
      fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
      fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
      fill3 <- sample(x = zero3,size = length(undetec3), replace = F)
      data_pboot <- data_p 
      data_pboot[fill1,1] <- undetec1
      data_pboot[fill2,2] <- undetec2
      data_pboot[fill3,3] <- undetec3
      bootx1<- rmultinom (n = 1, size = n[1], prob = data_pboot[,1])
      bootx2<- rmultinom (n = 1, size = n[2], prob = data_pboot[,2])
      bootx3<- rmultinom (n = 1, size = n[3], prob = data_pboot[,3])
      data_b<-cbind(bootx1,bootx2,bootx3)
      data_b <- data_b[rowSums(data_b)>0,] 
      colnames(data_b)<-c(colnames(data))
      data_b
      
    })
    
  }
  return(data_boot)
}

Chat1_f0Fun <-function(f1, f2, n) {
  if (f2 > 0) {
    f0 <- (n - 1) / n * f1^2 / (2 * f2)
    C <- 1 - f1 / n * (n-1)*f1/((n-1)*f1+2*f2)
    
  } else if (f2 == 0 & f1 != 0) {
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
    C <- 1 - f1 / n * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    
  } else {
    C <- 1
  }
  f0 <- ceiling(f0)
  return(c(C, f0))
}




miNEXT3 <- function(data, knots = 15, m1 = NULL,nboot = 0){
  if(is.null(colnames(data))){colnames(data) <- paste0("site",1:ncol(data))}
  n <- colSums(data)
  D <- apply(data,2, function(x) sum(x>0))
  
  #compute the mixture compositions
  m1 <- seq(0,n[1],length.out = knots) %>% round(.,0)
  m2 <- ((n[1]-m1)*n[2]/(n[2]+n[3])) %>% round(.,0)
  m3 <- ((n[1]-m1)*n[3]/(n[2]+n[3])) %>% round(.,0)
  
  #compute diveristy of each single community
  Each <- lapply(1:ncol(data),function(i){
    m <-  seq(0,max(n[i],n[1]),length.out = knots) %>% round(0)
    output <- iNEXT(x = data[,i],q = c(0,1,2),datatype = "abundance", size = m,se = F)$iNextEst
    output$qD[is.na(output$qD)] <- 0
    output %>% as_tibble() %>% select(m, q = order, Diversity=qD,method) %>% mutate(Community = names(n)[i])
  }) %>% do.call(rbind,.)
  Each$method[Each$method=="interpolated"] <- "rarefaction"
  Each$method[Each$method=="extrapolated"] <- "extrapolation"
  
  #compute the mixture diversity
  mix <- mix3(data,ms = cbind(m1,m2,m3)) %>% as.numeric()
  if(nboot>1){
    boot_sample<-Abun_CreatBootstrapSample(data=data,nboots = nboot)
    boot_se <- sapply(1:nboot,function(k){
      bdata <- boot_sample[[k]]
      bdata <- bdata[rowSums(bdata)>0,]
      mix3(bdata,ms = cbind(m1,m2,m3)) %>% as.numeric()
    }) %>% apply(., 1, sd)
  }else{
    boot_se <- rep(0,length(mix))
  }
  Mixture <- tibble(m1 = rep(m1,3), m2 = rep(m2,3), m3 = rep(m3,3),q = rep(c(0,1,2),each = length(m1)/3),
         Diversity = mix) %>% mutate(LCL = Diversity-qnorm(0.975)*boot_se,UCL = Diversity+qnorm(0.975)*boot_se,
                                     method = ifelse((m2<=n[2] & m3<=n[3]),"rarefaction","extrapolation"))
  return(list(Each = Each, Mixture = Mixture))
}

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
    q0 <- D0_rare_yhc(data[,1], data[,2],data[,3], m[1], m[2],m[3])
    q1 <- D_share_yhc(data[,1], data[,2],data[,3], m[1], m[2],m[3],1)
    q2 <- D2_yhc(data,m)
  }else if(sum(m[-1]>n[-1])==1){
    if(m[2] > n[2]){
      data <- data[,c(1,3,2)]
      datap <- datap[,c(1,3,2)]
      n <- n[c(1,3,2)]
      m <- m[c(1,3,2)]
    }
    q0 <- D_q0_ext1_yhc(data,datap,m,n)
    q1 <- D_q1_ext1_yhc(data,datap,m,n)
    q2 <- D2_yhc(data,m)
  }else{
    q0 <- D_q0_ext2_yhc(data,datap,m,n)
    q1 <- D_q1_ext2_yhc(data,datap,m,n)
    q2 <- D2_yhc(data,m)
  }
  return(c(q0,q1,q2))
}


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
#' @param m a vector of length 3 specifying number of individuals for each column of data 
#' @param q  a vector of diversity order, usually c(0,1,2)
#' @return a vector of mixture diversity of q = 0,1,2
#' @examples
#' data(Abudata)
#' data1<-Abudata$data
#' mm1<-Abudata$mm
#' Dq_in_yhc(data1[,1],data1[,2],data1[,3],mm1,0)
#' @export
Dq_in_yhc<-function(data,m,q){
  q1pls <- q[!c(q==0|q==2)]
  out1pls <- D_share_yhc(data[,1], data[,2],data[,3], m[1], m[2],m[3],q1pls)
  if(length(q1pls)>=2) out1pls <- t(out1pls)
  out0 <- D0_rare_yhc(data[,1], data[,2],data[,3], m[1], m[2],m[3])
  out2 <- D2_yhc(data = data,m)
  c(out0,out1pls,out2)
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

  tmp1 <- sum(mm*(mm-1)*pk/(sum(mm)^2))
  tmp2 <- 2*sapply(1:nrow(ch2),function(k){
    mm[ch2[k,]] %>% prod(.) * pkj[k]}
    ) %>% sum(.)/(sum(mm)^2)
  return((1/sum(mm) + tmp1 + tmp2)^(-1))

  # apply(mm, 1, subfun)
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
# q = 0, for two rarefaction and one extrapolation
D_q0_ext1_yhc<-function(data,datap,m,n){
  q0_tmp1 <- D0_rare_yhc(data[,1],data[,2],data[,3], m[1], m[2], n[3])
  h0_3_1hat <- h0_3_1hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2], m[3]-n[3],n[1],n[2],n[3])
  
  D_q0_ext1_value<-q0_tmp1+h0_3_1hat
  
  return(D_q0_ext1_value) 
}
# q = 0, for one rarefaction and two extrapolation
D_q0_ext2_yhc<-function(data,datap,m,n){
  q0_tmp1<-D0_rare_yhc(data[,1],data[,2],data[,3], m[1], n[2], n[3])
  h0_3_2hat <- h0_3_2hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2]-n[2], m[3]-n[3],n[1],n[2],n[3])
  D_q0_ext2_value<-q0_tmp1+h0_3_2hat
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
# q = 1, for two rarefaction and one extrapolation
D_q1_ext1_yhc<-function(data,datap,m,n){
  q1_tmp1 <- D_share_yhc(data[,1],data[,2],data[,3], m[1], m[2], n[3],1)
  h1_3_1hat <- h1_3_1hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2], m[3],n[1],n[2],n[3])
  
  D_q1_ext1_value<-q1_tmp1*exp(h1_3_1hat)
  return(D_q1_ext1_value)
}
# q = 1, for one rarefaction and two extrapolation
D_q1_ext2_yhc<-function(data,datap,m,n){
  q1_tmp1 <- D_share_yhc(data[,1],data[,2],data[,3], m[1], n[2], n[3],1)
  h1_3_2hat <- h1_3_2hat_cpp(datap[,1],datap[,2],datap[,3], m[1], m[2], m[3],n[1],n[2],n[3])
  D_q1_ext2_value<-q1_tmp1*exp(h1_3_2hat)
  return(D_q1_ext2_value)
}


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

ggmiNEXT3 <- function(x){
  nms <- unique(x$Each$Community)
  x$Mixture <- x$Mixture %>% mutate(Community = "Mixture") %>% select(m1,q,Diversity,LCL,UCL,method,Community)
  if(length(unique(x$Mixture$method))==2){
    border <- x$Mixture %>% filter(method=="rarefaction") %>% filter(m1 == min(m1))
    border$method <- "extrapolation"
    x$Mixture <- rbind(x$Mixture,border)
  }
  x$Each <- x$Each %>% mutate(LCL=Diversity,UCL=Diversity) %>% select(m1=m,q,Diversity,LCL,UCL,method,Community)
  x <- do.call(rbind,x)
  # output <- output %>% mutate(m1new = ifelse(Community == nms[2]|Community == nms[3],sum(data[,1])-m1, m1))
  x_h <- x %>% filter(method=="observed" & Community==nms[1]) %>% select(q,Diversity) %>% 
    group_by(q)
  x_p <- x %>% filter(method=="observed") 
  x$method[x$method=="observed"] <- "rarefaction" 
  x$method <- factor(x$method,levels = unique(x$method))
  pp = ggplot()+facet_grid(q~.,scales = "free_y")+
    geom_hline(data = x_h, aes(yintercept = Diversity), col = "darkgray", linetype = 3, size = 1.25)+
    geom_line(data = x, aes(x = m1, y = Diversity, col = Community, size = Community, lty = method))+
    geom_point(data = x_p, aes(x = m1, y = Diversity, col = Community),size = 3)+
    geom_ribbon(data = x, aes(x = m1, ymin = LCL, ymax = UCL, fill = Community),alpha = 0.4)+
    scale_size_manual(breaks=c(nms[1],nms[2],nms[3],"Mixture"), values=c(1, 1, 1, 1.3),
                      labels = c(nms[1],nms[2],nms[3],"Mixture"))+
    scale_color_manual(values = c("black", "#00BA38","#619CFF","#F8766D"),
                       breaks = c(nms[1],nms[2],nms[3],"Mixture"),
                       labels = c(nms[1],nms[2],nms[3],"Mixture"))+xlab("m1")+
    scale_fill_manual(values = c("black", "#00BA38","#619CFF","#F8766D"),
                      breaks = c(nms[1],nms[2],nms[3],"Mixture"))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(legend.position="bottom",legend.title=element_blank(),
          legend.text=element_text(size=10),legend.key.width  = unit(1.5,"cm"))+
    theme(plot.title = element_text(size=10, face="bold.italic",hjust = 0))+
    theme(axis.text.x= element_text(size = 10,colour = "black",margin=unit(c(0.2,0.2,0.5,0.5), "cm")),
          axis.text.y= element_text(size = 10,colour = "black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
          axis.title.x=element_text(size = 10,lineheight = 1.1),axis.title.y=element_text(size = 12),
          axis.ticks = element_line(size = 1.5,colour = "black"),axis.ticks.length = unit(-0.25,"line"))
  return(pp)
}



#' Mixture diversity of 3 communities
#' \code{miNEXT} computes the mixture diversity of 3 communities.
#' @param data a matrix/data.frame.
#' @param knots an integer specifying the number of points of m1. Default to be 15.
#' @param m1 a vector specifying the values of m1 where the mixture diversity will be computed.
#' @param nboot an integer specifying the number of bootstrap times to build confidence interval. 
#' Use 0 to skip bootstrap.
#' @return a list of 2 components: \code{$Each} a table of diversity of single community; \code{$Mixture} a 
#' table of mixture diversity.
#' @importFrom iNEXT iNEXT
#' @import dplyr 
#' @examples
#' data(Abudata)
#' data1 <- Abudata$data
#' result <- miNEXT3(data1)
#' @export
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

#' Plot function for output of miNEXT3.
#' \code{ggmiNEXT3} plot the output of miNEXT3 based on ggplot.
#' @param x output of miNEXT3.
#' @return a ggplot object
#' @import dplyr ggplot2
#' @examples
#' data(Abudata)
#' data1 <- Abudata$data
#' result <- miNEXT3(data1)
#' ggmiNEXT3(result)
#' @export
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




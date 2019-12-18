#######################################################################################################################
#' R code for 3-habitat mixture of 3-habitat mixture of rarefaction/extrapolation curves
#' @param data a vector/matrix/list of species abundance frequency
#' @param boot 
#' @param q  a numerical value for Dq
#' @return  a vector/matrix/list with species relative abundance or detection probability distribution
#' @examples
#' data(Abudata)
#' data1<-Abudata$data
#' mm1<-Abudata$mm
#' Dq_in(data1[,1],data1[,2],data1[,3],mm1,0)
#' @export

miNEXT3_Abu<-function(data=NULL,knots=10,size=NULL,nboots=0){
  result_abun<-NULL
  result_abun_CI<-NULL
  result_abun<-Abun3(data1=data,knots=knots,size=size)
  if(nboots>0){
    print("bootstrap start")
  }
}

Abun_ <- function(data1, knots = 10, size = NULL){
  x1 <- data1[, 1]
  x2 <- data1[, 2]
  x3 <- data1[, 3]
  n <- apply(data1,2,sum)
  D <- apply(data1,2, function(x) sum(x>0))
  
  m1 <- seq(0,n[1],length.out = 5) %>% round(.,0)
  m2 <- ((n[1]-m1)*n[2]/(n[2]+n[3])) %>% round(.,0)
  m3 <- ((n[1]-m1)*n[3]/(n[2]+n[3])) %>% round(.,0)
  
  #Compute each assemblage q=0,1,2
  each <- iNEXT(x = data1,q = c(0,1,2),datatype = "abundance",size = m1,se = F)$iNextEst %>% 
    sapply(., function(x) x$qD)
  each[is.na(each)] <- 0
  
  
  #Mixture
  # mix <- sapply(c(0,1,2), function(i){
  #   dq <-Dq_in(x1,x2,x3,cbind(m1,m2,m3),i)
  #   dq
  # })
  mm<-cbind(m1,m2,m3)
  mix_yhc <- Dq_in_yhc(x1,x2,x3,mm,c(0,1,2))
  
  
  #Mixture
  
  
}


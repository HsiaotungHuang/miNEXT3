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
#' @importFrom stats qnorm sd
#' @importFrom stats sd
#' @import dplyr 
#' @examples
#' \dontrun{
#' data(Abudata)
#' data1 <- Abudata$data
#' result <- miNEXT3(data1)
#' }
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
#' \dontrun{
#' data(Abudata)
#' data1 <- Abudata$data
#' result <- miNEXT3(data1)
#' ggmiNEXT3(result)
#' }
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
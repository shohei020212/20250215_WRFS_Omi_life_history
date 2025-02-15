
# Write by Noda on 2025.02.15

# 1. setup ---------------------------------------------------------------------

# package
{
  library("dplyr")
  require(ggplot2)
  require(flexmix)
}

# reset environment
rm(list=ls(all.names=T))

# 2. prepare raw data ----------------------------------------------------------

# AIC and BIC list
aic_bic <- data.frame(clu = c(),
                      aic = c(),
                      bic = c())

# number of k
kk <- 6


# 3. run flexmix ---------------------------------------------------------------

{
  # output name
  png_name <- paste("./output/figure_k", kk, ".png", sep = "")
  csv_name <- paste("./output/k", kk, "_data.csv", sep = "")
  
  # read data
  data <- read.csv("output/rawdata_filtered.csv")
  
  # filtering
  # remove individuals only captured once
  data <- data %>%
  # filter(id %in% names(which(table(id) > 1)),
  #        days_from_hatching < 1460) %>% 
    filter(treatment == 'early' | treatment == 'late' | treatment == 'control')
  
  ## von bertalanffy growth model ----------------------------------------------
  sbgm <- function(lmax,lmin,a1,theta,t0=0,t){
    (lt <- lmax-(lmax-lmin)*exp(-a1*(t-t0)-180*a1/pi*(cos(pi*(t0+theta)/180)-cos(pi*(t+theta)/180))))
  }
  
  # functions
  llsbgm <- function(para,data){
    lmax <- exp(para[1])
    lmin <- exp(para[2])
    a1 <- exp(para[3])
    theta <- para[4]
    t <- data[,1]
    lobs <- data[,2]
    ll <- -log(dnorm(log(lobs),log(sbgm(lmax,lmin,a1,theta,t=t)),sd=sd(log(lobs))))
    sum(ll)
  }
  
  llsbgm2 <- function(para,data){
    lmax <- exp(para[1])
    lmin <- exp(para[2])
    a1 <- exp(para[3])
    theta <- para[4]
    t <- data[,1]
    lobs <- data[,2]
    ll <- -log(dnorm(log(lobs),log(sbgm(lmax,lmin,a1,theta,t=t)),sd=sd(log(lobs))))
    ll
  }
  
  sbgmplot <- function(par,t,ylim,col="black",main="",xlab="",ylab="",type="l"){
    lmax <- exp(par[1])
    lmin <- exp(par[2])
    a1 <- exp(par[3])
    theta <- par[4]
    plot(t,sbgm(lmax,lmin,a1,theta,t=t),ylim=ylim,col=col,main=main,xlab=xlab,ylab=ylab,type=type,
         font.lab=2,cex.main=2.5,cex.lab=2.4,cex.axis=2,mgp=c(4.5, 1, 0), lwd=4, las=1)
  }
  
  ipara <- c(log(max(data$fl)),log(min(data$fl)),log(1.7E-4),0)
  
  lmax <- exp(ipara[1])
  lmin <- exp(ipara[2])
  a1 <- exp(ipara[3])
  theta <- ipara[4]
  
  t <- seq(0,max(data$days_from_hatching),by=1)
  
  llsbgm(ipara,data[,c("days_from_hatching","fl")])
  summary(llsbgm2(ipara,data[,c("days_from_hatching","fl")]))
  
  # single cluster (Models that do not assume growth classes)
  res <- optim(ipara,llsbgm,data=data[,c("days_from_hatching","fl")])
  
  #AIC
  sum(res$value)+2*5
  #BIC
  sum(res$value)+5*log(nrow(data))
  
  lmax <- exp(res$par[1])
  lmin <- exp(res$par[2])
  a1 <- exp(res$par[3])
  theta <- res$par[4]
  
  png('./output/figure_growth_Curve.png', width = 1000, height = 600)
  par(family="Arial",mar=c(7, 8, 5, 5),mgp=c(4.5, 1, 0))
  plot(t,sbgm(lmax,lmin,a1,theta,t=t),type="l",lwd = 3, col=8,ylim=c(0,max(data$fl)), 
       font.lab=2,
       xlab="No. of days since January-1",ylab="Fork length (mm)",
       cex.lab=3,cex.axis=2,las=1)
  points(fl~days_from_hatching,data=data,pch=16,cex=1,col=rgb(0, 0, 0, alpha=0.1))
  dev.off()
  
  data$predfl <- sbgm(lmax,lmin,a1,theta,t=data$days_from_hatching)
  data$resid <- log(data$fl)-log(data$predfl)
  
  # model
  m1 <- FLXMRglm(resid~1,family="gaussian")
  
  set.seed(1)
  
  mp2 <- flexmix(.~.|id, data=data, model=list(m1), k=kk, 
                 control=list(iter=200,minprior=0,nrep=3,classify="random"))
  
  # show summary
  summary(mp2)
  
  data$cluster <- clusters(mp2)
  
  k <- sort(unique(data$cluster))
  
  ipara <- res$par
  uid <- sort(unique(data$id))
  
  ll <- as.numeric()
  oll <- 1
  set.seed(1)
  while(abs(sum(ll)-sum(oll))>0.1){
    # M step
    oll <- ll
    ll <- as.numeric()
    estpara <- as.numeric()
    for(cl in k){
      res <- optim(ipara,llsbgm,data=data[data$cluster==cl,c("days_from_hatching","fl")])
      estpara <- rbind(estpara,res$par)
      data[,paste("ll",cl,sep="-")]<- llsbgm2(res$par,data[,c("days_from_hatching","fl")])
      ll <- c(ll,res$value)
    }
  }
  data$ocluster <- data$cluster
  # E step
  for(cid in uid){
    dmy <- data[data$id==cid,]
    cll <- apply(dmy[,paste("ll",k,sep="-")],2,sum)
    data[data$id==cid,]$cluster <- which(cll==min(cll))
  }
  print(c(sum(ll),sum(oll)))
  
  
  data$ocluster <- data$cluster
  
  #AIC
  aic <- sum(ll)+2*max(k)*5
  #BIC
  bic <- sum(ll)+max(k)*5*log(nrow(data))
  
  # make new data frame
  add_data <- data.frame(clu = kk,
                         aic = aic,
                         bic = bic)
  aic_bic <- rbind(aic_bic, add_data)
  
  # color list
  cb = brewer.pal(kk, "RdYlBu")

  # sorting
  par(mfrow=c(1,max(k)))
  
  # make new data frame
  ocluster <- c()
  max_len <- c()
  clu_max_len <- data.frame(ocluster,max_len)
  
  for(ocluster in k){
    tmp = data[data$ocluster==ocluster,]
    tmp = tmp[tmp$days_from_hatching<300,]
    max_len = max(tmp$fl)
    # make new data frame
    add_data <- data.frame(ocluster, max_len)
    clu_max_len <- rbind(clu_max_len, add_data)
  }
  

  # sort by max length
  clu_max_len <- clu_max_len[order(clu_max_len$max_len, decreasing=T),]
  # create a sequential list
  x <- 1:kk
  # make new data frame
  add_data <- data.frame(x)
  names(add_data) <- c("ocluster_rename")
  # bind
  clu_max_len <- cbind(clu_max_len, add_data)
  # merge
  data <- left_join(data, clu_max_len, by=c("ocluster"="ocluster"))
  
  # Change the order of the growth curves
  # merge a sequential list
  estpara2 <- cbind(estpara, x)
  # make new data frame
  estpara2 <- data.frame(estpara2)
  # merge
  estpara2 <- left_join(estpara2, clu_max_len, by=c("x"="ocluster"))
  # sort
  estpara2 <- estpara2[order(estpara2$ocluster_rename),]
  # remove columns
  estpara2 <- estpara2[,-5]
  # convert matrix
  estpara2 <- matrix(as.matrix(estpara2), nrow(estpara2), ncol(estpara2))
  
  
  ### plot growth history
  kkk <- max(k)
  if (max(k) == kk) {
    num_lis <- seq(1,kk,1)
  } else {
    num_lis <- seq(1,kkk,1)
  }
  # plot
  png(png_name, width = 1200, height = 900, pointsize = 15, bg = "white", res=100)
  par(family="Arial",mar=c(7, 7, 5, 5))
  for(i in num_lis) {
    p1 <- sbgmplot(estpara2[i,],t,c(0,max(data$fl)),col=cb[i],xlab="Days from hatching",ylab="Fork length (mm)") +
      points(fl~days_from_hatching,data=data[data$ocluster_rename==i,],bg=cb[i],pch=21,cex=0.8,lwd=0.5)
    grid(col = "gray", lty = 1)
    legend("bottomright",
           title="growth cluster",
           cex=1.5,
           legend = k,
           fill = cb,
           horiz=TRUE)
    par(new=T)
  }
  graphics.off()
  
  # out to csv
  data2 <- aggregate(ocluster_rename~id,data,unique)
  write.csv(data2,csv_name, row.names = FALSE)
}

# ...----


# 4. check BIC -----------------------------------------------------------------

# plot 
g <- ggplot(aic_bic,aes(x=clu,y=bic)) +
  geom_line(linewidth = 2) +
  geom_point(size = 5) + 
  xlab("Number of cluster") +
  ylab("BIC") +
  theme_bw() +
  theme(strip.text=element_text(size=16, family="sans", face="bold"),
        axis.title=element_text(size=16, family="sans", face="bold"),
        axis.text=element_text(size=16, family="sans"),
        strip.background = element_blank(),
        panel.background = element_rect("white"))
g

# save png format
ggsave(g,
       file = paste("output/figure_flexmix_bic.png", sep = ""),
       width = 4,
       height = 3)

# save csv
write.csv(aic_bic, "output/flexmix_bic.csv")



# end --------------------------------------------------------------------------

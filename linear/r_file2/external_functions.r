xprob=as.matrix(cbind(expand.grid(0:3,0:2,0:2,0:2),1,1,1,1))
x1p=c(0.2,0.36,0.34,0.1);
x2p=c(0.28,0.55,0.17);
x3p= c(0.7408182,0.2222455, 1-0.7408182-0.2222455)
x4p= c(0.8187308,0.1637462,1-0.8187308-0.1637462)
for(i in 1:108){
  xprob[i,5]=x1p[xprob[i,1]+1]*x2p[xprob[i,2]+1]*x3p[xprob[i,3]+1]*x4p[xprob[i,4]+1]
  
}
xprob<-xprob[order(xprob[,1],xprob[,2],xprob[,3],xprob[,4]),]

tmpl_pe<-function(betaxin1,cut_q){
  betaxin=betaxin1
  xprob[,6]=betaxin[1]+as.matrix(xprob[,1:4])%*%(betaxin[-1]) #phix
  xprob=xprob[order(xprob[,6]),] # order by phix
  xprob[,7]=cumsum(xprob[,5]) # cumsum
  for(i in 2:108){ # group
    if(xprob[i,7]<cut_q[1]|xprob[i-1,7]<cut_q[1]){
      xprob[i,8]=1
    }
    else if(xprob[i,7]<cut_q[2]|xprob[i-1,7]<cut_q[2]){
      xprob[i,8]=2
    }
    else if(xprob[i,7]<cut_q[3]|xprob[i-1,7]<cut_q[3]){
      xprob[i,8]=3
    }
    else{ xprob[i,8]=4}
    
  }
  tmp_l=xprob[order(xprob[,1],xprob[,2],xprob[,3],xprob[,4]),8]
  # pe
  pe_test=c()
  for(j in 1:4){
    dat_temp=xprob[which(xprob[,8]==j),]
    pe_test=c(pe_test, sum(dat_temp[,6]*dat_temp[,5])/sum(dat_temp[,5]))
    
  }
  return(list(tmp_l,pe_test))
}



#xprob[order(xprob[,1],xprob[,2],xprob[,3],xprob[,4]),8]
#pe

external_infor<-function(sig_y,target_beta,NN=50000,m=1000, var_if){
  
  
  #External Data
  
  betaxin=matrix(nrow=m,ncol = 5)
  betaxin_var=matrix(nrow=m,ncol = 5*5)
  
  tmplist=matrix(nrow=m,ncol = 108)
  pe=matrix(nrow=m,ncol = 4)
  tmplist2=matrix(nrow=m,ncol = 108)
  pe2=matrix(nrow=m,ncol = 4)
  
  ref_sam_freq=matrix(0,nrow=m,ncol = 108)

  cutX1 <- c(-Inf, qnorm(c(0.2,0.56,0.9),0,1), Inf)
  cutX2 <- c(-Inf, qnorm(c(0.28,0.83),0,1), Inf)
  
  for(i in 1:m){
    
    X1 <- rnorm(NN,0,1)
    X2 <- rnorm(NN,0,1)
    X3 <- rpois(NN,0.3)
    X4 <- rpois(NN,0.2)
    
    X1 <- as.numeric(cut(X1, breaks=cutX1, right=T)) - 1
    X2 <- as.numeric(cut(X2, breaks=cutX2, right=T)) - 1
    X3[which(X3 > 2)] <- 2
    X4[which(X4 > 2)] <- 2
    
    Q <- cbind(1,X1,X2,X3,X4)
    theta <- c(-2,0.1,-0.1,0.1,0.1)
    sig_logZ <- 0.6
    u_logZ <- Q%*%theta 
    logZ <- rtnorm(NN, mean=u_logZ, sd=sig_logZ, lower=-Inf, upper=0)#rnorm(1,u_Z,sig_Z)
    Z <- exp(logZ)
    Zcat <- as.numeric(cut(Z, breaks=c(-100,c(1:9)/10,1.1),right=F)) - 1
    
    
    alpha <-target_beta[1] ; beta_X <- target_beta[2:5]; beta_Z <-target_beta[6]
    tmp <- (cbind(1,X1,X2,X3,X4,Zcat)%*%c(alpha,beta_X,beta_Z))
    Y <- tmp+rnorm(NN,sd=sig_y)
    
    Dat.external <- data.frame(Y=Y, X1=X1,X2=X2,X3=X3,X4=X4,Z=Zcat)
    
    
    
    #####Phi(X)####################################################
    lm.external <- lm(Y ~ X1+X2+X3+X4,data=Dat.external)
    
    betaxin_t <- as.vector(lm.external$coefficients)
    betaxin[i,]<-betaxin_t
    if(var_if==T){
      betaxin_var[i,]<-as.vector(vcov(lm.external))
    }
    
    
    tmplist[i,]<-tmpl_pe(betaxin_t,cut_q = c(0.25,0.5,0.75))[[1]]
    pe[i,]<-tmpl_pe(betaxin_t,cut_q = c(0.25,0.5,0.75))[[2]]
    
    tmplist2[i,]<-tmpl_pe(betaxin_t,cut_q = c(0.5,0.7,0.9))[[1]]
    pe2[i,]<-tmpl_pe(betaxin_t,cut_q = c(0.5,0.7,0.9))[[2]]
    
    sampleref= cbind(X1,X2,X3,X4)
    ref_sam_freq1<-table(sampleref%*%c(27,9,3,1)+1)/NN
    ref_ind<-as.numeric(row.names(ref_sam_freq1))
    
    ref_sam_freq[i,ref_ind]<-ref_sam_freq1
    
  }
  

  return(list(betaxin,tmplist,pe,betaxin_var,tmplist2,pe2,ref_sam_freq))
  
}

sample_ref<-function(nr){
  
  cutX1 <- c(-Inf, qnorm(c(0.2,0.56,0.9),0,1), Inf)
  cutX2 <- c(-Inf, qnorm(c(0.28,0.83),0,1), Inf)
  
  X1 <- rnorm(nr,0,1)
  X2 <- rnorm(nr,0,1)
  X3 <- rpois(nr,0.3)
  X4 <- rpois(nr,0.2)
  
  X1 <- as.numeric(cut(X1, breaks=cutX1, right=T)) - 1
  X2 <- as.numeric(cut(X2, breaks=cutX2, right=T)) - 1
  X3[which(X3 > 2)] <- 2
  X4[which(X4 > 2)] <- 2
  
  return( cbind(X1,X2,X3,X4))
  
  
}

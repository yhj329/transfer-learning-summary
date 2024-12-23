simu_glm2<-function(type_ini,target_beta =c(-2.4,-0.5,0.4,0.3,0.65,0.1),random_target,source_beta,
                   ns=2000,N_r=1000,nt,xtol=1.0e-4,al_type,vartype=0,nr,nr_type){
  # vartype: 0 no additional term; 1 using diag from target; 
  # 2 using full matrix from target; 3 using asy formula
  # random_target==T random external otherwise fixed
  # type_ini   initial_value  1:(betaxin 0)  2:(betaxin 0.2) otherwise mle
  # nr reference sample  nr_type 0:before repetition 1: in repetition
  #------------------------------------------------------------------
  
  # results store
  # estimate
  beta_mle<-matrix(ncol=6+6,nrow=N_r)
  beta_my<-matrix(ncol=6,nrow=N_r)
  beta_cmle<-matrix(ncol=6,nrow=N_r)
  beta_cmle2<-matrix(ncol=6,nrow=N_r)
  
  asyd_mle<-matrix(ncol=6+6,nrow=N_r)
  asyd_cml<-matrix(ncol=6,nrow=N_r)
  asyd_cml2<-matrix(ncol=6,nrow=N_r)
  asyd_my<-matrix(ncol=6,nrow=N_r)
  actnn=matrix(0,nrow=N_r,ncol=8)
  actnn2=matrix(0,nrow=N_r,ncol=8)
  
  # source data
  cutX1 <- c(-Inf, qnorm(c(0.05,0.21,0.77),0,1), Inf)
  cutX2 <- c(-Inf, qnorm(c(0.23,0.82),0,1), Inf)
  sampleref=sample_ref(nr)
  for(i in 1:N_r){
    
    X1 <- rnorm(ns,0,1)
    X2 <- rnorm(ns,0,1)
    X3 <- rpois(ns,0.5)
    X4 <- rpois(ns,0.2)
    
    X1 <- as.numeric(cut(X1, breaks=cutX1, right=T)) - 1
    X2 <- as.numeric(cut(X2, breaks=cutX2, right=T)) - 1
    X3[which(X3 > 2)] <- 2
    X4[which(X4 > 2)] <- 2
    
    Q <- cbind(1,X1,X2,X3,X4)
    theta <- c(-2,0.1,-0.1,0.1,0.1)
    sig_logZ <- 0.6
    u_logZ <- Q%*%theta 
    logZ <- rtnorm(ns, mean=u_logZ, sd=sig_logZ, lower=-Inf, upper=0)#rnorm(1,u_Z,sig_Z)
    Z <- exp(logZ)
    Zcat <- as.numeric(cut(Z, breaks=c(-100,c(1:9)/10,1.1),right=F)) - 1
    
    alpha <-source_beta[1] ; beta_X <-source_beta[2:5]; beta_Z <-source_beta[6]
    tmp <- (cbind(1,X1,X2,X3,X4,Zcat)%*%c(alpha,beta_X,beta_Z))
    Y <- rbinom(ns,size=1,prob=1/(1+exp(-tmp)))
    
    
    # mle and gammahat
    res_mle<-glm(Y~X1+X2+X3+X4+Zcat,family = binomial)
    beta_mle[i,1:6]<-res_mle$coefficients
    
    
    gammahat_tot<-nlm(f =function(p)gammah_likelihood(p,cbind(1,X1,X2,X3,X4),Z),p = c(0,0,0,0,0,1),hessian = T)
    gammahat<-gammahat_tot$estimate # gamma estimate
    beta_mle[i,7:12]<-gammahat
    asyv_gamma=solve(gammahat_tot$hessian) #asy var gamma
    asyd_mle[i,7:12]<-sqrt(diag(asyv_gamma))  #asy se gamma
    
    
    
    # cmle function(data,gammahat,pe,tmpl) data is y x zc
    
    
    
    # target matrix
    xprob2=matrix(0,nrow=nrow(xprob),ncol=4+1+10)
    xprob2[,1:5]<-xprob[,1:5]
    xprob2[,6]<-0
    for(j in 0:9){
      
      # zcgiven X
      pvec=(pnorm(log((j+1)/10),mean=gammahat[1]+xprob2[,1:4]%*%gammahat[2:5],sd=gammahat[6])-
              pnorm(log(j/10),mean=gammahat[1]+xprob2[,1:4]%*%gammahat[2:5],sd=gammahat[6]))/pnorm(0,mean=gammahat[1]+xprob2[,1:4]%*%gammahat[2:5],sd=gammahat[6])
      xprob2[,j+5+1]<-pvec # conditional distribution of z given x 
      
      #beta[1]+xprob[,1:4]%*%beta[2:5]+beta[6]*j*pvec
    }
    
    if(random_target==T){
      betaxin_1<-betaxin_total[[1]][i,]
      pe_cml=betaxin_total[[3]][i,]
      tmpl_cml=betaxin_total[[2]][i,]
      pe_cml2=betaxin_total[[6]][i,]
      tmpl_cml2=betaxin_total[[5]][i,]
      varxin=matrix(betaxin_total[[4]][i,],nrow=5)
    }
    else{ pe_cml=betaxin_total[[3]][1,]
    tmpl_cml=betaxin_total[[2]][1,]
    pe_cml2=betaxin_total[[6]][1,]
    tmpl_cml2=betaxin_total[[5]][1,]
    betaxin_1<-betaxin_total[[1]][1,]
    varxin=matrix(betaxin_total[[4]][1,],nrow=5)
    }
    
    # initial
    initial_value=res_mle$coefficients
    if(type_ini==1){
      initial_value=c(betaxin_1,0)
    }else if(type_ini==2){
      initial_value=c(betaxin_1,0.2)
    }
    
    
    
    
    
    # umat<-diag(rep(1,6))
    
    
    fisher_ma2<-fisher_info(betahat =res_mle$coefficients,datax = cbind(1,X1,X2,X3,X4),datazcat = Zcat)
    asy_mle_beta<-solve(fisher_ma2)/ns
    asyd_mle[i,1:6]<-sqrt(diag(asy_mle_beta))
    
    
    
    
    
    #-----------------my------------------------------------------------
    betazt=as.numeric(res_mle$coefficients[6])
    
    # sample reference sample
    if(nr_type==1){
      sampleref=sample_ref(nr)
    }
    
    ref_sam_freq<-table(sampleref%*%c(27,9,3,1)+1)/nr
    ref_ind<-as.numeric(row.names(ref_sam_freq))
    
    xprob3<-xprob2
    xprob3[,5]<-0
    xprob3[ref_ind,5]<-ref_sam_freq
    
    #betamy<-rootSolve::multiroot(f=function(beta)
    # est_f2(beta,xprobdata = xprob2,xref =sampleref ,betaxin=betaxin_1,betaz=betazt),start = betaxin_1)$root
    betamy<-rootSolve::multiroot(f=function(beta)est_f(beta,xprobdata = xprob2,betaxin=betaxin_1,betaz=betazt),start = betaxin_1)$root
    
    
    
    
    
    beta_my[i,1:5]<-betamy
    beta_my[i,6]<-betazt
    
    #agamma
    Agamma=A_gamma(c(betamy,betazt),xprobdata = xprob3,gammah = gammahat)
    
    #ax
    Ax=rootSolve::gradient(function(beta)est_f(beta,xprobdata = xprob3,betaxin=betaxin_1,betaz=betazt),x = betamy)
    
    #az
    Az= rootSolve::gradient(function(betaz)est_f(beta=betamy,xprobdata = xprob3,betaxin=betaxin_1,betaz=betaz),x = betazt)
    
    # axin
    Axin=rootSolve::gradient(function(betaxin)est_f(beta=betamy,xprobdata = xprob3,betaxin,betaz=betazt),x = betaxin_1)
    # vartype: 0 no additional term; 1 using diag from target; 
    # 2 using full matrix from target; 3 using asy formula
    term_xin=0
    if(vartype==0){
      term_xin=0
    }else if(vartype==1){
      vxin=diag(diag(varxin))
      term_xin=(Axin)%*%vxin%*%(Axin)/nt
    }else if(vartype==2){
      vxin=varxin
      term_xin=(Axin)%*%vxin%*%(Axin)/nt
    }else if(vartype==3){
      term_xin=var_h(beta=betamy,xprobdata = xprob3,betaxin = betaxin_1,betaz = betazt)/nt
    }
    
    
    
    
    asyv_my<-solve(Ax)%*%(Az%*%t(Az)*solve(fisher_ma2)[6,6]/ns +Agamma%*%asyv_gamma%*%t(Agamma)+term_xin)%*% t(solve(Ax) )
    asyd_my[i,1:5]<-sqrt( diag(asyv_my  )   )
    asyd_my[i,6]<-sqrt(solve(fisher_ma2)[6,6]/ns)
    # Bgamma<-B_gamma(xprobdata =xprob2,gammahat,betazt)
    # asyd_my[i,]<-c( sqrt(diag( Bgamma%*%asyv_gamma%*%t(Bgamma)+exz%*%t(exz)*asy_mle_beta[6,6] )), sqrt(asy_mle_beta[6,6]))#
  }  
  
  return(list(beta_mle=beta_mle,beta_cmle=beta_cmle,beta_cmle2=beta_cmle2,beta_my=beta_my,asyd_mle=asyd_mle,asyd_cml=asyd_cml,asyd_cml2=asyd_cml2,asyd_my=asyd_my,actnn=actnn,actnn2=actnn2))
  
  
}
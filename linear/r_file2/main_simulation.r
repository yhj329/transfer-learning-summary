
simu_linear<-function(betaxin_total,sig_y =0.8,type_ini,target_beta =c(-2.4,-0.5,0.4,0.3,0.65,0.5),random_target=T,source_beta,ns=2000,nt,N_r=1000,xtol=1.0e-4,vartype=0,al_type,nr){
  # vartype: 0 no additional term;  3 using asy formula
  # random_target==T random external otherwise fixed
  # type_ini   initial_value  1:(betaxin 0)  2:(betaxin 0.2) otherwise mle
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
  
  actnn<-0
  
  # source data
  cutX1 <- c(-Inf, qnorm(c(0.05,0.21,0.77),0,1), Inf)
  cutX2 <- c(-Inf, qnorm(c(0.23,0.82),0,1), Inf)

  
  
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
    Y <- tmp+rnorm(ns,sd=sig_y)
    
    # mle and gammahat
    res_mle<-lm(Y~X1+X2+X3+X4+Zcat)
    beta_mle[i,1:6]<-res_mle$coefficients
    
    
    gammahat_tot<-nlm(f =function(p)gammah_likelihood(p,cbind(1,X1,X2,X3,X4),Z),p = c(0,0,0,0,0,1),hessian = T)
    gammahat<-gammahat_tot$estimate # gamma estimate
    beta_mle[i,7:12]<-gammahat
    asyv_gamma=solve(gammahat_tot$hessian) #asy var gamma
    asyd_mle[i,7:12]<-sqrt(diag(asyv_gamma))  #asy se gamma
    
    
    sigma_mle=sqrt(mean(res_mle$residuals^2))
    
    # cmle function(data,gammahat,pe,tmpl) data is y x zc
    
    # target matrix
    xprob2=xprob
    xprob2[,6]<-0
    vz<-rep(0,nrow(xprob))
    for(j in 0:9){
      
      # zcgiven X
      pvec=(pnorm(log((j+1)/10),mean=gammahat[1]+xprob2[,1:4]%*%gammahat[2:5],sd=gammahat[6])-
              pnorm(log(j/10),mean=gammahat[1]+xprob2[,1:4]%*%gammahat[2:5],sd=gammahat[6]))/pnorm(0,mean=gammahat[1]+xprob2[,1:4]%*%gammahat[2:5],sd=gammahat[6])
      xprob2[,6]<-xprob2[,6]+(j*pvec) # conditional mean of z given x 
      vz<-vz+((j^2)*pvec)
      #beta[1]+xprob[,1:4]%*%beta[2:5]+beta[6]*j*pvec
    }
    
    Ezz=vz%*%xprob2[,5]
    
    
    
    
    
    if(random_target==T){
      betaxin_1<-betaxin_total[[1]][i,]
      pe_cml=betaxin_total[[3]][i,]
      tmpl_cml=betaxin_total[[2]][i,]
      pe_cml2=betaxin_total[[6]][i,]
      tmpl_cml2=betaxin_total[[5]][i,]
      ref_sam_freq=betaxin_total[[7]][i,]
      
      # varxin=matrix(betaxin_total[[4]][i,],nrow=5)
    }
    else{ pe_cml=betaxin_total[[3]][1,]
    tmpl_cml=betaxin_total[[2]][1,]
    betaxin_1<-betaxin_total[[1]][1,]
    pe_cml2=betaxin_total[[6]][1,]
    tmpl_cml2=betaxin_total[[5]][1,]
    ref_sam_freq=betaxin_total[[7]][1,]
    #varxin=matrix(betaxin_total[[4]][1,],nrow=5)
    }
    
    # initial
    initial_value=res_mle$coefficients
    if(type_ini==1){
      initial_value=c(betaxin_1,0)
    }else if(type_ini==2){
      initial_value=c(betaxin_1,0.2)
    }
    
    
    
    #------- estimate cml1----------------------------
    res_cml <- nloptr(x0=initial_value,
                      eval_f=function(beta)likelihood_beta(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat),
                      eval_grad_f=function(beta)likelihood_beta_score(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat),
                      eval_g_ineq = function(beta)constr_cml(beta,xprobdata=xprob2,pe=pe_cml,tmpl =tmpl_cml),
                      eval_jac_g_ineq = function(beta)Cderivative2(betahat=beta,xprobdata =xprob2,tmpl = tmpl_cml ),
                      opts = list("algorithm"=al_type,"xtol_rel"=xtol)
    )
    beta_cmle[i,]<-res_cml$solution         # "NLOPT_LN_COBYLA"   "NLOPT_LD_MMA" "NLOPT_LD_SLSQP"
    
    
    # asyd
    
    
    cder<-Cderivative(betahat =res_cml$solution,xprobdata = xprob2,pe=pe_cml,tmpl =tmpl_cml)
    if(length(cder)==1){
      umat<-diag(rep(1,6))
      actnn=actnn+1
    }
    else{umat<-Null(t(cder))}
    
    
    fisher_ma<-fisher_info(cbind(1,X1,X2,X3,X4),Zcat,sigma_mle)
    vmat=umat%*%solve(t(umat)%*%fisher_ma%*%umat)%*%t(umat)
    Wma=cbind(1,X1,X2,X3,X4,Zcat)
    middle_I<-matrix(0,nrow=6,ncol=6)
    for(k in 1:ns){
      constt=as.numeric((Y[k]-t(beta_cmle[i,])%*%Wma[k,])^2)
      middle_I=middle_I+Wma[k,]%*%t(Wma[k,])*constt/ns
      
    }
    middle_I<-middle_I/(sigma_mle^4)
    asy_cml=vmat%*%middle_I%*%vmat/ns
    
    asyd_cml[i,]<-sqrt(diag(asy_cml))
    asy_mle_beta<-solve(fisher_ma)/ns
    asyd_mle[i,1:6]<-sqrt(diag(asy_mle_beta))
    
    
    #------- estimate cml2----------------------------
    res_cml2 <- nloptr(x0=initial_value,
                      eval_f=function(beta)likelihood_beta(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat),
                      eval_grad_f=function(beta)likelihood_beta_score(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat),
                      eval_g_ineq = function(beta)constr_cml(beta,xprobdata=xprob2,pe=pe_cml2,tmpl =tmpl_cml2),
                      eval_jac_g_ineq = function(beta)Cderivative2(betahat=beta,xprobdata =xprob2,tmpl = tmpl_cml2 ),
                      opts = list("algorithm"=al_type,"xtol_rel"=xtol)
    )
    beta_cmle2[i,]<-res_cml2$solution         # "NLOPT_LN_COBYLA"   "NLOPT_LD_MMA" "NLOPT_LD_SLSQP"
    
    
    # asyd
    
    
    cder<-Cderivative(betahat =res_cml2$solution,xprobdata = xprob2,pe=pe_cml2,tmpl =tmpl_cml2)
    if(length(cder)==1){
      umat<-diag(rep(1,6))
      actnn=actnn+1
    }
    else{umat<-Null(t(cder))}
    

    vmat=umat%*%solve(t(umat)%*%fisher_ma%*%umat)%*%t(umat)
    Wma=cbind(1,X1,X2,X3,X4,Zcat)
    middle_I<-matrix(0,nrow=6,ncol=6)
    for(k in 1:ns){
      constt=as.numeric((Y[k]-t(beta_cmle2[i,])%*%Wma[k,])^2)
      middle_I=middle_I+Wma[k,]%*%t(Wma[k,])*constt/ns
      
    }
    middle_I<-middle_I/(sigma_mle^4)
    asy_cml2=vmat%*%middle_I%*%vmat/ns
    
    asyd_cml2[i,]<-sqrt(diag(asy_cml2))

    
    #------------------------ my-------------------------------
    
    
    xprob3<-xprob
    xprob3[,5]<-ref_sam_freq
    
    Exx=vxx(xprob3)
    
    betazt=as.numeric(res_mle$coefficients[6])
    exz=as.vector(t(xprob2[,6]*xprob3[,5])%*%cbind(1,xprob2[,1:4]))
    betamyx<-betaxin_1-solve(Exx)%*%(exz)*betazt
    beta_my[i,1:5]<-betamyx
    beta_my[i,6]<-betazt
    
    
    if(vartype==0){
      v_xin=0
    }else if(vartype==3){
      v_xin=sigma_mle^2*solve(Exx)
      residuals_temp=(betamyx[1]-betaxin_1[1]+xprob[,1:4]%*%(betamyx[2:5]-betaxin_1[2:5]))
      residuals_xin=residuals_temp^2+2*residuals_temp*xprob2[,6]*betazt+vz*betazt^2
      
      E_xin_middle=t(cbind(1,xprob2[,1:4]))%*%( ( (as.vector(xprob2[,5])*as.vector(residuals_xin) )
                                                  *cbind(1,xprob2[,1:4])) )
      v_xin=v_xin+solve(Exx)%*%E_xin_middle %*%solve(Exx)
    }
    
    
    
    
    Bgamma<-B_gamma(xprobdata =xprob3,gammahat,betazt)
    asyd_my[i,]<-c( sqrt(diag(solve(Exx)%*%  Bgamma%*%asyv_gamma%*%t(Bgamma)%*%solve(Exx)+solve(Exx)%*%exz%*%t(exz)%*%solve(Exx)*asy_mle_beta[6,6]+v_xin/nt )), sqrt(asy_mle_beta[6,6]))#
  }  
  
  return(list(beta_mle=beta_mle,beta_cmle=beta_cmle,beta_cmle2=beta_cmle2,beta_my=beta_my,asyd_mle=asyd_mle,asyd_cml=asyd_cml,asyd_my=asyd_my,asyd_cml2=asyd_cml2))
  
  
}

gammah_likelihood<-function(gammah,datax,dataz){
  return(-sum(log(dtnorm(x=log(dataz),mean=datax%*%gammah[1:5],sd=gammah[6],lower=-Inf,upper=0))))
  
}

#(glm)
likelihood_beta<-function(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat){
  
  tmp=1/(1+exp(-datax%*%beta[1:5]-datazc*beta[6])) 
  return( -mean(datay*log(tmp)+(1-datay)*log(1-tmp)) )   # negative likelihood
}

likelihood_beta_score<-function(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat){
  tmp=exp(-datax%*%beta[1:5]-datazc*beta[6]) 
  
  return(-t(cbind(datax,datazc))%*%((datay*tmp/(1+tmp))-(1-datay)/(1+tmp))/nrow(datax)) # negative likelihood
}


constr_cml<-function(beta,xprobdata,pe=betaxin_total[[3]][1,],tmpl=betaxin_total[[2]][1,]){
  pre_pe=c()
  meany_x<-rep(0,nrow(xprobdata))
  for(j in 1:nrow(xprobdata)){
    tmp=exp(-beta[1]-as.numeric(xprobdata[j,1:4]%*%beta[2:5])-c(0:9)*beta[6] )
    meany_x[j]<- sum(xprobdata[j,6:15]/(1+tmp))
  }
  
  for(k in 1:4){
    index_pe=which(tmpl==k)
    pre_pe=c(pre_pe,sum(meany_x[index_pe]*xprobdata[index_pe,5])
             /sum(xprobdata[index_pe,5]))
    
  }
  return(c(pre_pe-pe-0.1*abs(pe),-pre_pe+pe-0.1*abs(pe) ) )
  
}

Cderivative<-function(tol_con=1/100000,xprobdata,betahat,pe,tmpl){
  act_pe=which(abs(constr_cml(betahat,xprobdata,pe,tmpl))<tol_con)
  if(length(act_pe)==0){
    return(1)}
  c_der=Cderivative2(betahat =betahat,xprobdata = xprobdata,tmpl = tmpl )[act_pe,]
  return(c_der)
}

Cderivative2<-function(betahat,xprobdata,tmpl){
  
  meanyd_x<-matrix(0,nrow=nrow(xprobdata),ncol=length(betahat))
  for(j in 1:nrow(xprobdata)){
    tmp=exp(-betahat[1]-as.numeric(xprobdata[j,1:4]%*%betahat[2:5])-c(0:9)*betahat[6] )
    meanyd_x[j,1:5]<- sum(xprobdata[j,6:15]*tmp/(1+tmp)^2)*c(1,xprobdata[j,1:4])
    meanyd_x[j,6]<- sum(xprobdata[j,6:15]*(0:9)*tmp/(1+tmp)^2)
  }
  
  c_der=matrix(0,nrow=8,ncol=6)
  for(i in 1:4){
    
    # vec=rep(0,6)
    index=which(tmpl==i)
    # for(j in index){
    #   vec=vec+xprobdata[j,5]*c(1,xprobdata[j,1:4],xprobdata[j,6])
    # }
    c_der[i,]= t(xprobdata[index,5])%*%meanyd_x[index,]/sum(xprobdata[index,5])
    
  }
  c_der[5:8,]=-c_der[1:4,]
  return(c_der)
}


# fisher information
fisher_info<-function(betahat,datax,datazcat){
  dataxz=cbind(datax,datazcat)
  tmp=exp(-dataxz%*%betahat)
  tmp=as.vector(tmp/(1+tmp)^2)
  # broadcast?
  cov_inf=t(tmp*dataxz)%*%dataxz/(nrow(dataxz))
  return(cov_inf)
  
}

Middle_I<-function(datay,betahat,datax,datazcat){
  dataxz=cbind(datax,datazcat)
  tmp=exp(-dataxz%*%betahat)
  tmp=as.vector( ((datay*tmp)/(1+tmp)-(1-datay)/(1+tmp))^2)
  # broadcast?
  middleI=t(tmp*dataxz)%*%dataxz/(nrow(dataxz))
  return(middleI)
  
}

est_f<-function(beta,xprobdata,betaxin,betaz){
  meanyp_x<-rep(0,nrow(xprobdata))
  for(j in 1:nrow(xprobdata)){
    tmp=exp(-beta[1]-as.numeric(xprobdata[j,1:4]%*%beta[2:5])-c(0:9)*betaz)
    
    meanyp_x[j]<- sum(xprobdata[j,6:15]/(1+tmp))#
  }
  tmpxin=exp(-betaxin[1]-(xprobdata[,1:4]%*%betaxin[2:5]))
  return(t(cbind(1,xprobdata[,1:4]))%*%( (meanyp_x-1/(1+tmpxin))*xprobdata[,5]))
}




est_f2<-function(beta,xprobdata,xref,betaxin,betaz){
  total_sum=rep(0,5)
  for(i in 1:nrow(xref)){
    tmp=exp(-beta[1]-as.numeric(xref[i,1:4]%*%beta[2:5])-c(0:9)*betaz)
    ind=as.numeric(xref[i,1:4]%*%c(27,9,3,1)+1)
    meanyp_x=sum(xprobdata[ind,6:15]/(1+tmp))#
    
    tmpxin=exp(-betaxin[1]-(xref[i,1:4]%*%betaxin[2:5]))
    total_sum=total_sum+c(1,xref[i,1:4])*(meanyp_x-1/(1+tmpxin))
    
  }
  
  return(total_sum)
}







pgamma_zj<-function(gammah,zj,x){
  result=pnorm(log((j+1)/10),mean=x%*%gammah[1:5],sd=gammah[6])-pnorm(log(j/10),mean=x%*%gammah[1:5],sd=gammah[6])
  return(result/pnorm(0,mean=x%*%gammah[1:5],sd=gammah[6]))
  
}
der_pgamma_zj<-function(betahat,gammah,x){
  tmp=exp(-as.numeric(x%*%betahat[1:5])-(0:9)*betahat[6])
  dj=dnorm(log(c(0:10)/10),mean=x%*%gammah[1:5],sd=gammah[6])
  pj=pnorm(log(c(0:10)/10),mean=x%*%gammah[1:5],sd=gammah[6])
  alphaj=(log(c(0:10)/10)-as.numeric(x%*%gammah[1:5]))/gammah[6]
  ad_j=alphaj*dj
  ad_j[1]<-0
  const1=-( pj[11]*(dj[2:11]-dj[1:10]) -dj[11]*(pj[2:11]-pj[1:10]))/(pj[11]^2) 
  const1=sum(const1/(1+tmp))
  const2=-( pj[11]*(ad_j[2:11]-ad_j[1:10])-ad_j[11]*(pj[2:11]-pj[1:10]))/(pj[11]^2) 
  const2=sum(const2/(1+tmp))
  return(c(const1,const2))
  #return(list(dj,pj,alphaj,ad_j))
}
A_gamma<-function(betahat,xprobdata,gammah){
  bga<-matrix(0,nrow=5,ncol=6)
  for(i in 1:108){
    der_ga_x<-der_pgamma_zj(betahat,gammah,c(1,xprobdata[i,1:4]))
    bga=bga+xprobdata[i,5]*c(1,xprobdata[i,1:4])%*%t(c(der_ga_x[1]*c(1,xprobdata[i,1:4]),der_ga_x[2]))
  }
  return(bga)
}
#return(der_ga_x)

var_h<-function(beta,xprobdata,betaxin,betaz){
  meanyp_x<-rep(0,nrow(xprobdata))
  for(j in 1:nrow(xprobdata)){
    tmp=exp(-beta[1]-as.numeric(xprobdata[j,1:4]%*%beta[2:5])-c(0:9)*betaz)
    
    meanyp_x[j]<- sum(xprobdata[j,6:15]/(1+tmp))#
  }
  tmpxin=exp(-betaxin[1]-(xprobdata[,1:4]%*%betaxin[2:5]))
  
  scalar_v=as.vector(   (meanyp_x*(1-2/(1+tmpxin)) +1/(1+tmpxin)^2
  )*xprobdata[,5])
  
  return(t(cbind(1,xprobdata[,1:4]))%*% diag(scalar_v) %*%cbind(1,xprobdata[,1:4]) )
}





gammah_likelihood<-function(gammah,datax,dataz){
  return(-sum(log(dtnorm(x=log(dataz),mean=datax%*%gammah[1:5],sd=gammah[6],lower=-Inf,upper=0))))
  
}

#(ls)
likelihood_beta<-function(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat){
  return(sum((datay-datax%*%beta[1:5]-datazc*beta[6])^2)/nrow(datay)) # negative likelihood
}

likelihood_beta_score<-function(beta,datay=Y,datax=cbind(1,X1,X2,X3,X4),datazc=Zcat){
  return(-2*t(cbind(datax,datazc))%*%(datay-datax%*%beta[1:5]-datazc*beta[6])/nrow(datay)) # negative likelihood
}

constr_cml<-function(beta,xprobdata,pe=betaxin_total[[3]][1,],tmpl=betaxin_total[[2]][1,]){
  pre_pe=c()
  for(k in 1:4){
    index_pe=which(tmpl==k)
    pre_pe=c(pre_pe,sum((beta[1]+xprobdata[index_pe,1:4]%*%beta[2:5]+xprobdata[index_pe,6]*beta[6])*xprobdata[index_pe,5])
             /sum(xprobdata[index_pe,5]))
    
  }
  return(c(pre_pe-pe-0.1*abs(pe),-pre_pe+pe-0.1*abs(pe)))
  
}

Cderivative<-function(tol_con=1/100000,xprobdata,betahat,pe,tmpl){
  act_pe=which(abs(constr_cml(betahat,xprobdata,pe,tmpl))<tol_con)
  act_pe_n=length(act_pe)
  if(act_pe_n==0){
    return(1)
  }
  c_der=matrix(nrow=act_pe_n,ncol=6)
  for(i in 1:act_pe_n){
    actcol=act_pe[i]
    if(actcol>4){
      sign_pre=-1
      intv=actcol-4}
    else {sign_pre=1;intv=actcol}
    
    # vec=rep(0,6)
    index=which(tmpl==intv)
    # for(j in index){
    #   vec=vec+xprobdata[j,5]*c(1,xprobdata[j,1:4],xprobdata[j,6])
    # }
    vec=c(sum(xprobdata[index,5]*1),xprobdata[index,5]%*%xprobdata[index,1:4],
          sum(xprobdata[index,5]*xprobdata[index,6])  )
    
    c_der[i,]=sign_pre*vec/sum(xprobdata[index,5])
    
  }
  return(c_der)
}

Cderivative2<-function(betahat,xprobdata,tmpl){
  
  c_der=matrix(0,nrow=8,ncol=6)
  for(i in 1:4){
    
    # vec=rep(0,6)
    index=which(tmpl==i)
    # for(j in index){
    #   vec=vec+xprobdata[j,5]*c(1,xprobdata[j,1:4],xprobdata[j,6])
    # }
    vec=c(sum(xprobdata[index,5]*1),xprobdata[index,5]%*%xprobdata[index,1:4],
          sum(xprobdata[index,5]*xprobdata[index,6])  )
    
    c_der[i,]=vec/sum(xprobdata[index,5])
    
  }
  c_der[5:8,]=-c_der[1:4,]
  return(c_der)
}


# fisher information
fisher_info<-function(datax,datazcat,sigm){
  dataxz=cbind(datax,datazcat)
  cov_inf=t(dataxz)%*%dataxz/(nrow(dataxz)*sigm^2)
  return(cov_inf)
  
}

pgamma_zj<-function(gammah,zj,x){
  result=pnorm(log((j+1)/10),mean=x%*%gammah[1:5],sd=gammah[6])-pnorm(log(j/10),mean=x%*%gammah[1:5],sd=gammah[6])
  return(result/pnorm(0,mean=x%*%gammah[1:5],sd=gammah[6]))
  
}
der_pgamma_zj<-function(gammah,x){
  dj=dnorm(log(c(0:10)/10),mean=x%*%gammah[1:5],sd=gammah[6])
  pj=pnorm(log(c(0:10)/10),mean=x%*%gammah[1:5],sd=gammah[6])
  alphaj=(log(c(0:10)/10)-as.numeric(x%*%gammah[1:5]))/gammah[6]
  ad_j=alphaj*dj
  ad_j[1]<-0
  const1=-( pj[11]*(dj[2:11]-dj[1:10]) -dj[11]*(pj[2:11]-pj[1:10]))/(pj[11]^2) 
  const1=sum(const1*(0:9))
  const2=-( pj[11]*(ad_j[2:11]-ad_j[1:10])-ad_j[11]*(pj[2:11]-pj[1:10]))/(pj[11]^2) 
  const2=sum(const2*(0:9))
  return(c(const1,const2))
  #return(list(dj,pj,alphaj,ad_j))
}
B_gamma<-function(xprobdata,gammah,betazhat){
  bga<-matrix(0,nrow=5,ncol=6)
  for(i in 1:108){
    der_ga_x<-der_pgamma_zj(gammah,c(1,xprobdata[i,1:4]))
    bga=bga+xprobdata[i,5]*c(1,xprobdata[i,1:4])%*%t(c(der_ga_x[1]*c(1,xprobdata[i,1:4]),der_ga_x[2]))
  }
  return(-bga*betazhat)
  #return(der_ga_x)
}
vxx<-function(xprobdata){
  
  vxx=t(cbind(1,xprobdata[,1:4]))%*%(as.vector(xprobdata[,5])*cbind(1,xprobdata[,1:4]))
  
  return(vxx)
}


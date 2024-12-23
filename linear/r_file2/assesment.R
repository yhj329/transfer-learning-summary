

assesment<-function(result,targetbeta){
  
  total_comp<-data.frame(c("S1","","","","",""),
                         c("$\\beta_{T,X_0}$","$\\beta_{T,X_1}$","$\\beta_{T,X_2}$","$\\beta_{T,X_3}$","$\\beta_{T,X_4}$","$\\beta_{T,Z^c}$"),                    targetbeta,
                         t(t(colMeans(result$beta_mle[,1:6]))-targetbeta),
                         colMeans(result$asyd_mle[,1:6]),
                         sqrt(diag(var(result$beta_mle[,1:6]))),
                         colMeans(abs(t(t(result$beta_mle[,1:6])-targetbeta))<qnorm(0.975)*result$asyd_mle[,1:6]),
                         t(t(colMeans(result$beta_cmle[,1:6]))-targetbeta),
                         colMeans(result$asyd_cml[,1:6]),
                         sqrt(diag(var(result$beta_cmle[,1:6]))),
                         colMeans(abs(t(t(result$beta_cmle[,1:6])-targetbeta))<qnorm(0.975)*result$asyd_cml[,1:6]),
                         
                         
                         t(t(colMeans(result$beta_cmle2[,1:6]))-targetbeta),
                         colMeans(result$asyd_cml2[,1:6]),
                         sqrt(diag(var(result$beta_cmle2[,1:6]))),
                         colMeans(abs(t(t(result$beta_cmle2[,1:6])-targetbeta))<qnorm(0.975)*result$asyd_cml2[,1:6]),
                         
                         t(t(colMeans(result$beta_my[,1:6]))-targetbeta),
                         colMeans(result$asyd_my[,1:6]),
                         sqrt(diag(var(result$beta_my[,1:6]))),
                         colMeans(abs(t(t(result$beta_my[,1:6])-targetbeta))<qnorm(0.975)*result$asyd_my)
                         
  )
  rownames(total_comp)<-c("$\\beta_{T,X_0}$","$\\beta_{T,X_1}$","$\\beta_{T,X_2}$","$\\beta_{T,X_3}$","$\\beta_{T,X_4}$","$\\beta_{T,Z^c}$")
  
  return(total_comp)
  
}
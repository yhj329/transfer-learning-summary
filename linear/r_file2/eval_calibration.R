# 该函数用于评估calibration图像 输入真实模型f_r 预测模型 f_p（x，z） quantile 个数r_q 随机产生N full target样本(x,y,z)

function(f_p,r_q,f_r,N){
  cutX1 <- c(-Inf, qnorm(c(0.2,0.56,0.9),0,1), Inf)
  cutX2 <- c(-Inf, qnorm(c(0.28,0.83),0,1), Inf)
  
  X1 <- rnorm(N,0,1)
  X2 <- rnorm(N,0,1)
  X3 <- rpois(N,0.3)
  X4 <- rpois(N,0.2)
  
  X1 <- as.numeric(cut(X1, breaks=cutX1, right=T)) - 1
  X2 <- as.numeric(cut(X2, breaks=cutX2, right=T)) - 1
  X3[which(X3 > 2)] <- 2
  X4[which(X4 > 2)] <- 2
  
  Q <- cbind(1,X1,X2,X3,X4)
  theta <- c(-2,0.1,-0.1,0.1,0.1)
  sig_logZ <- 0.6
  u_logZ <- Q%*%theta 
  logZ <- rtnorm(N, mean=u_logZ, sd=sig_logZ, lower=-Inf, upper=0)#rnorm(1,u_Z,sig_Z)
  Z <- exp(logZ)
  Zcat <- as.numeric(cut(Z, breaks=c(-100,c(1:9)/10,1.1),right=F)) - 1
  
  predict_p<-f_p(Q,Zcat)
  cut_p<-quantile(predict_p,r_q)
  
  #category
  mean_int<-c(1:r_q)
  mean_y_int<-c(1:r_q)
  label_p<-cut(predict_p, breaks=cut_p, right=T)
  
  for(i in 1:r_q){
    index_int<-which(label_p==i)
    mean_int[i]=mean(predict_p[index_int])
    mean_y_int[i]=mean(y[index_int])
  }
  plot(mean_y_int,mean_int)
  lines()
  
}



# 接下来
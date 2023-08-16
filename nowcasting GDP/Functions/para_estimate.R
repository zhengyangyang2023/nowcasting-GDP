EM_DFM_SS_block_idioQARMA_restrMQ<-function(X,Par){
  
  #X=x
  
  thresh = 1e-4
  r = Par$r
  p = Par$p
  max_iter = Par$max_iter
  i_idio = Par$i_idio
  R_mat = Par$Rconstr #4行5列
  q = Par$q
  nQ = Par$nQ
  blocks = Par$blocks
  
  ### Prepara??o dos dados
  
  TT <- dim(X)[1] #397
  N <- dim(X)[2] #25
  
  ### Standardise X
  Mx = colMeans(X,na.rm=T)
  Wx = sapply(1:N,function(x) sd(X[,x],na.rm = T))
  xNaN <- (X-kronecker(t(Mx),rep(1,TT)))/kronecker(t(Wx),rep(1,TT)) 
  
  ### Initial conditions
  
  # Removing missing values (for initial estimator)
  optNaN<-list()
  optNaN$method = 2; # Remove leading and closing zeros
  optNaN$k = 3;
  
  res_InitCond<-InitCond(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio);
  
  A<-res_InitCond$A #53阶
  C<-res_InitCond$C #25行53列
  Q<-res_InitCond$Q #53阶
  R<-res_InitCond$R #25阶
  Z_0<-res_InitCond$initZ #53*1，全为0
  V_0<-res_InitCond$initV #53*53
  
  # some auxiliary variables for the iterations
  previous_loglik = -Inf
  num_iter = 0
  LL = -Inf
  converged = F
  
  # y for the estimation is WITH missing data
  y = t(xNaN) #xNaN为397*25
  
  #--------------------------------------------------------------------------
  #THE EM LOOP
  #--------------------------------------------------------------------------
  
  #The model can be written as
  #y = C*Z + e;
  #Z = A*Z(-1) + v
  #where y is NxT, Z is (pr)xT, etc
  
  #remove the leading and ending nans for the estimation
  optNaN$method = 3
  y_est <- remNaNs_spline(xNaN,optNaN)
  y_est_indNaN<-t(y_est$indNaN) #25*385的缺失情况矩阵
  y_est<-t(y_est$X) #25*385
  
  
  num_iter = 0
  loglik_aux = -Inf
  
  # 使用EM算法更新A, C, Q, R, Z_0, V_0
  while ((num_iter < max_iter) & !converged){
    
    
    res_EMstep = EMstep(y_est, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks)
    # res_EMstep <- list(C_new, R_new, A_new, Q_new, Z_0, V_0, loglik)
    
    C = res_EMstep$C_new;
    R = res_EMstep$R_new;
    A = res_EMstep$A_new;
    Q = res_EMstep$Q_new;
    Z_0<-res_EMstep$Z_0
    V_0<-res_EMstep$V_0
    loglik<-res_EMstep$loglik
    
    # updating the user on what is going on
    # if((num_iter%% 50)==0&&num_iter>0){message(num_iter,"th iteration: \nThe loglikelihood went from ",round(loglik_aux,4)," to ",round(loglik,4))}
    loglik_aux <- loglik
    
    # Checking convergence
    if (num_iter>2){
      res_em_converged = em_converged(loglik, previous_loglik, thresh,1)
      converged<-res_em_converged$converged
      decreasse<-res_em_converged$decrease
    }
    
    LL <- cbind(LL, loglik)
    previous_loglik <- loglik
    num_iter <-  num_iter + 1
  }
  
  # final run of the Kalman filter
  res_runKF = runKF(y, A, C, Q, R, Z_0, V_0)
  Zsmooth<-t(res_runKF$xsmooth) #398*53。原数据后面整行缺失的已经预测好。无缺失值。
  x_sm <- Zsmooth[2:dim(Zsmooth)[1],]%*%t(C) #397行25列
  # x_sm不就是下面的FF乘以C吗
  
  Res<-list()
  Res$x_sm <- x_sm #数据的平滑值(标准化后)397行25列
  Res$X_sm <- kronecker(t(Wx),rep(1,TT))*x_sm + kronecker(t(Mx),rep(1,TT))#数据的平滑值(未标准化)397行25列
  Res$FF <- Zsmooth[2:dim(Zsmooth)[1],]#397行53列，状态变量的平滑值。
  
  #--------------------------------------------------------------------------
  #  Loading the structure with the results
  #--------------------------------------------------------------------------
  Res$C <- C;
  Res$R <- R;
  Res$A <- A;
  Res$Q <- Q;
  Res$Mx <- Mx;#X的每列的均值(去除na)
  Res$Wx <- Wx;#X的每列的标准差(去除na)
  Res$Z_0 <- Z_0;
  Res$V_0 <- V_0;
  Res$r <- r;
  Res$p <- p;
  Res$loglik <- loglik;
  
  return(Res)
  
  
}
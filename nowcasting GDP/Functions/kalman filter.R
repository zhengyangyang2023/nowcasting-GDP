FIS <- function(Y,Z,R,TT,Q,S){
  
  m<-dim(S$Am)[1]
  nobs<-dim(S$Am)[2]
  
  S$AmT           = matlab::zeros(m,nobs+1)
  S$PmT           = array(0,c(m,m,nobs+1))
  S$AmT[,nobs+1] <- S$AmU[,nobs+1]
  S$PmT[,,nobs+1] <- S$PmU[,,nobs+1]
  
  S$PmT_1<-array(0,c(m,m,nobs))
  S$PmT_1[,,nobs] <- (diag(m)-S$KZ)%*%TT%*%S$PmU[,,nobs]
  
  if(length(S$Pm[,,nobs])==1){pinv = 1/S$Pm[,,nobs]}else{pinv<-corpcor::pseudoinverse(S$Pm[,,nobs])}
  
  J_2 <- S$PmU[,,nobs]%*%t(TT)%*%pinv
  
  for (t in nobs:1){ 
    PmU <- S$PmU[,,t]
    Pm1 <- S$Pm[,,t]
    P_T <- S$PmT[,,t+1]
    P_T1 <- S$PmT_1[,,t]
    
    J_1 <- J_2
    
    S$AmT[,t] <- S$AmU[,t] + J_1%*%(S$AmT[,t+1] - TT%*%S$AmU[,t])
    S$PmT[,,t] <- PmU + J_1%*%(P_T - Pm1)%*%t(J_1) 
    
    if(t>1){
      if(length(S$Pm[,,t-1])==1){pinv = 1/S$Pm[,,t-1]}else{pinv<-corpcor::pseudoinverse(S$Pm[,,t-1])}
      J_2 <- S$PmU[,,t-1]%*%t(TT)%*%pinv
      S$PmT_1[,,t-1] = PmU%*%t(J_2)+J_1%*%(P_T1-TT%*%PmU)%*%t(J_2)
    }
  }
  
  return(S)
  
}

SKF <-function(Y,Z,R,TT,Q,A_0,P_0){
  
  n <- dim(Z)[1]
  m <- dim(Z)[2]
  nobs  <- dim(Y)[2]
  
  S<-list()
  S$Am <- array(NA,c(m,nobs)) #m行nobs列的NA值
  S$Pm <- array(NA,c(m,m,nobs))
  S$AmU <- array(NA,c(m,nobs+1))
  S$PmU <- array(NA,c(m,m,nobs+1))
  S$loglik <- 0
  
  Au <- A_0;  # A_0|0;
  Pu <- P_0;  # P_0|0
  
  S$AmU[,1]    = Au;
  S$PmU[,,1]  = Pu;
  
  
  
  for(t in 1:nobs){
    
    A <- TT%*%Au;
    P <- TT%*%Pu%*%t(TT) + Q;
    P <-  0.5*(P+t(P))
    
    # handling the missing data
    res_MissData <- MissData(y=Y[,t],C=Z,R=R)
    
    y_t <- res_MissData$y
    Z_t <- res_MissData$C
    R_t <- res_MissData$R
    L_t <- res_MissData$L
    
    # if(is.null(y_t)){
    if(sum(is.na(y_t))==length(y_t)){
      Au <- A
      Pu <- P
    } else {
      
      
      if(!is.matrix(Z_t)){
        Z_t<-t(as.matrix(Z_t)) #又加上了转置
      }
      
      
      PZ  <- P%*%t(Z_t)
      iF  <- solve(Z_t%*%PZ + R_t)
      PZF <- PZ%*%iF
      
      V <- y_t - Z_t%*%A
      Au <- A  + PZF%*%V
      Pu <- P  - PZF%*%t(PZ)
      Pu <-  0.5*(Pu+t(Pu))
      S$loglik <- S$loglik + 0.5*(log(det(iF))  - t(V)%*%iF%*%V) #对数似然函数越大越好
    }
    
    S$Am[,t] <- A
    S$Pm[,,t] <- P
    
    # Au = A_t|t   & Pu = P_t|t
    
    S$AmU[,t+1] <- Au
    S$PmU[,,t+1] <- Pu
  } # t
  
  if(sum(is.na(y_t))==length(y_t)){
    S$KZ <- matlab::zeros(m,m)
  }else{
    S$KZ <- PZF%*%Z_t
  }
  
  return(S)
  
}

MissData <- function(y,C,R){
  
  ix <- !is.na(y) # y为向量
  e  <- diag(length(y))
  L  <- e[,ix]
  
  y <-    y[ix]
  if(dim(C)[1]==1){C = C[ix]}else{C  <-  C[ix,]}
  R  <-  R[ix,ix]
  
  return(list(y=y,C=C,R=R,L=L))
  
}

runKF <- function(y, A, C, Q, R, x_0, Sig_0){
  
  # x_0 = Z_0; Sig_0 = V_0
  S <- SKF(Y=y,Z=C,R=R,TT=A,Q=Q, A_0=x_0, P_0=Sig_0);
  S <- FIS(Y=y,Z=C,R=R,TT=A,Q=Q,S=S);
  
  xsmooth <- S$AmT;
  Vsmooth <- S$PmT;
  VVsmooth <- S$PmT_1;
  loglik <- S$loglik;
  
  return(list(xsmooth = xsmooth,Vsmooth = Vsmooth,VVsmooth = VVsmooth,loglik = loglik))  
  
}
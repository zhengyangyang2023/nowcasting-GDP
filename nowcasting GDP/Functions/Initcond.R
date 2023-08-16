InitCond<-function(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio){
  
  x<-xNaN
  Rcon<-R_mat
  
  # library(magic)
  
  pC = dim(Rcon)[2] #5
  ppC = max(p,pC) #5
  n_b = dim(blocks)[2] #4
  
  res_remNaNs_spline <- remNaNs_spline(x,optNaN) 
  xBal<- res_remNaNs_spline$X 
  indNaN <- res_remNaNs_spline$indNaN 
  TT <- dim(xBal)[1] # time T of series
  N <- dim(xBal)[2] # number of series
  NM <- N-nQ # number of monthy frequency series：23
  
  xNaN = xBal
  for(i in 1:N){
    xNaN[indNaN[,i],i] <- NA
  }  
  
  # Initialize model coefficient output
  C = {}
  A = {}
  Q = {}
  initV = {}
  
  res = xBal #填充后的数据
  resNaN = xNaN #未填充的数据
  
  # Set the first observations as NaNs: For quarterly-monthly aggreg. scheme
  indNaN[1:pC-1,] <- T #pC=5，将缺失情况中的前4行判为缺失
  
  #r为4个1, 4为block个数
  for(i in 1:n_b){    # loop for each block
    #i=1
    r_i<-r[i]
    
    ########################
    # Observation equation #
    ########################
    
    C_i = matlab::zeros(N,r_i*ppC) # Initialize state variable matrix helper.25行5列
    idx_i = matlab::find(blocks[,i]) # returns the series loaded in block i
    # Note that the variables have been reshuffled so as to have all quarterly at the last columns of X
    idx_iM = idx_i[idx_i<NM+1];   # index for monthly variables
    idx_iQ = idx_i[idx_i>NM];     # index for quarterly variables
    
    
    eig<-eigen(cov(as.matrix(res[,idx_iM],ncol=length(idx_iM))))
    v<-eig$vectors[,1:r_i]
    d<-eig$values[1:r_i] 
    C_i[idx_iM,1:r_i] = v 
    f = as.matrix(res[,idx_iM])%*%as.matrix(v)
    
    
    for(kk in 0:(max(p+1,pC)-1)){ 
      if(kk == 0){
        FF<-f[(pC-kk):(dim(f)[1]-kk),]
      }else{
        FF <- cbind(FF,f[(pC-kk):(dim(f)[1]-kk),]) 
      } 
    }
    
    
    
    Rcon_i = kronecker(Rcon,diag(r_i)) 
    q_i = kronecker(q,matlab::zeros(r_i,1)) 
    ff = FF[,1:(r_i*pC)] 
    
    for(j in idx_iQ){     
      # j=1
      xx_j = resNaN[pC:dim(resNaN)[1],j]
      if(sum(!is.na(xx_j)) < dim(ff)[2]+2){ 
        xx_j = res[pC:dim(res)[1],j]
      }
      ff_j = ff[!is.na(xx_j),] 
      xx_j = xx_j[!is.na(xx_j)]
      iff_j = solve(t(ff_j)%*%ff_j) 
      Cc = iff_j%*%t(ff_j)%*%xx_j 
      Cc = Cc - iff_j%*%t(Rcon_i)%*%solve(Rcon_i%*%iff_j%*%t(Rcon_i))%*%(Rcon_i%*%Cc-q_i);#Rcon_i是约束矩阵。这是有约束的最小二乘回归系数
      C_i[j,1:(pC*r_i)] <- t(Cc)
    }
    ff = rbind(matlab::zeros(pC-1,pC*r_i),ff) 
    res = res - ff%*%t(C_i) 
    resNaN = res
    for(i_aux in 1:dim(indNaN)[2]){
      resNaN[indNaN[,i_aux],i_aux] <- NA 
    }
    C <- cbind(C,C_i) 
    
    #######################    
    # Transition Equation #
    #######################  
    
    z <- FF[,1:r_i] 
    Z <- FF[,(r_i+1):(r_i*(p+1))] 
    A_i = t(matlab::zeros(r_i*ppC,r_i*ppC)) 
    A_temp = solve(t(Z)%*%Z)%*%t(Z)%*%z 
    A_i[1:r_i,1:(r_i*p)] <- t(A_temp)
    A_i[(r_i+1):dim(A_i)[1],1:(r_i*(ppC-1))] <- diag(r_i*(ppC-1))
    
    Q_i = matlab::zeros(ppC*r_i,ppC*r_i) 
    e = z  - Z%*%A_temp         # VAR residuals
    Q_i[1:r_i,1:r_i] = cov(e);  # VAR covariance matrix
    
    initV_i = matrix(solve(diag((r_i*ppC)^2)-kronecker(A_i,A_i))%*%c(Q_i) ,r_i*ppC,r_i*ppC);
    
    
    
    #A Q initV准对角阵堆叠#
    if(is.null(A)){
      A<-A_i
    }else{
      A <- magic::adiag(A,A_i)  #准对角阵
    }
    
    if(is.null(Q)){
      Q<-Q_i
    }else{
      Q <- magic::adiag(Q,Q_i)  
    }
    
    if(is.null(initV)){
      initV<-initV_i
    }else{
      initV <- magic::adiag(initV,initV_i)  
    }
  }
  
  R = diag(apply(resNaN, 2, stats::var, na.rm = T))
  
  eyeN = diag(N)
  eyeN<-eyeN[,i_idio] 
  
  # Initial conditions
  C <- cbind(C,eyeN) 
  
  ii_idio = which(i_idio) 
  n_idio = length(ii_idio) 
  B = matlab::zeros(n_idio)
  S = matlab::zeros(n_idio)
  
  BM <- matlab::zeros(n_idio)
  SM <- matlab::zeros(n_idio)
  
  # Loop for monthly variables
  for (i in 1:n_idio){
    # Set observation equation residual covariance matrix diagonal
    R[ii_idio[i],ii_idio[i]] <- 1e-04 
    
    # Subsetting series residuals for series i
    res_i <- resNaN[,ii_idio[i]] 
    res_i<-res_i[!is.na(res_i)] 
    
    # Linear regression: AR 1 process for monthly series residuals
    BM[i,i] = solve(t(res_i[1:(length(res_i)-1)])%*%res_i[1:(length(res_i)-1)])%*%t(res_i[1:(length(res_i)-1)])%*%res_i[2:length(res_i)] 
    SM[i,i] = stats::var(res_i[2:length(res_i)]-res_i[1:(length(res_i)-1)]*BM[i,i])
  }
  
  # blocks for covariance matrices
  initViM = diag(1/diag(diag(dim(BM)[1])-BM^2))%*%SM; #23*23。1-(BM的对角线元素平方) 的倒数构成对角阵,再乘以SM
  
  if(!nQ==0){ 
    C<-cbind(C,rbind(matlab::zeros(NM,3*nQ),t(kronecker(diag(nQ),c(1,1,1)))))
    Rdiag<-diag(R)
    sig_e <- Rdiag[(NM+1):N]/19 
    Rdiag[(NM+1):N] <- 1e-04
    R <- diag(Rdiag) # Covariance for obs matrix residuals
    
    # for BQ, SQ
    rho0 <- 0.1
    temp <- matlab::zeros(3)
    temp[1,1] <- 1
    
    # BQ and SQ
    BQ <- kronecker(diag(nQ),rbind(cbind(rho0,matlab::zeros(1,2)),cbind(diag(2),matlab::zeros(2,1)))) #直积，A%X%B,A的每一个元素乘到B上
    if(length(sig_e)>1){ 
      SQ = kronecker(diag((1-rho0^2)*sig_e),temp)
    }else{
      SQ = kronecker((1-rho0^2)*sig_e,temp)
    }
    
    initViQ = matrix(solve(diag((3*nQ)^2)-kronecker(BQ,BQ))%*%c(SQ),3*nQ,3*nQ)
    
  }else{
    BQ <- NULL
    SQ <- NULL
    initViQ <- NULL
  }  
  
  A1<-magic::adiag(A,BM,BQ) #准对角阵，对角块为A，BM，BQ
  Q1<-magic::adiag(Q, SM, SQ)
  
  A<-A1 # Observation matrix
  Q<-Q1 # Residual covariance matrix
  
  # Initial conditions
  initZ = matlab::zeros(dim(A)[1],1); # States
  initV = magic::adiag(initV, initViM, initViQ) # Covariance of states
  
  return(list(A = A, C = C, Q = Q, R = R, initZ = initZ, initV = initV))
}
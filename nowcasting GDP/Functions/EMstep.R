EMstep <- function(y = NULL, A = NULL, C = NULL, Q = NULL, R = NULL, Z_0 = NULL, V_0 = NULL, 
                   r = NULL, p = NULL, R_mat = NULL, q = NULL, nQ = NULL, i_idio = NULL, blocks = NULL){
  
  n <- dim(y)[1] #25
  TT <- dim(y)[2] #385
  nM <- n - nQ #23
  pC <- dim(R_mat)[2] #5
  ppC <- max(p,pC) #5
  n_b <- dim(blocks)[2] #4
  
  # Compute the (expected) sufficient statistics for a single Kalman filter sequence.
  
  #Running the Kalman filter with the current estimates of the parameters
  res_runKF = runKF(y, A, C, Q, R, Z_0, V_0);
  
  
  Zsmooth<-res_runKF$xsmooth #53行386列
  Vsmooth<-res_runKF$Vsmooth
  VVsmooth<-res_runKF$VVsmooth
  loglik<-res_runKF$loglik
  
  A_new <- A
  Q_new <- Q
  V_0_new <- V_0
  
  #更新A，Q，V0的前4个对角阵
  for(i in 1:n_b){
    
    r_i <- r[i]
    rp <- r_i*p
    if(i == 1){
      rp1 <- 0*ppC
    }else{
      rp1 <- sum(r[1:(i-1)])*ppC #i=4时，为3r*5=15r
    }
    A_i <- A[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
    Q_i <- Q[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
    
    if(rp==1){
      EZZ <- t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)] +
        sum(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),2:dim(Vsmooth)[3]])  # E(Z'Z)
      EZZ_BB <- t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) %*% Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)] + 
        sum(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),1:(dim(Vsmooth)[3]-1)]) #E(Z(-1)'Z_(-1))
      EZZ_FB <- t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)] + 
        sum(VVsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),]) #E(Z'Z_(-1))
      
    }else{
      EZZ <- (Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) +
        apply(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),2:dim(Vsmooth)[3]],c(1,2),sum)  # E(Z'Z)。
      EZZ_BB <- (Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) %*% t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) + 
        apply(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),1:(dim(Vsmooth)[3]-1)],c(1,2),sum) #E(Z(-1)'Z_(-1))
      EZZ_FB <- (Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) + 
        apply(VVsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),],c(1,2),sum) #E(Z'Z_(-1))
    }
    
    
    A_i[1:r_i,1:rp] <- EZZ_FB[1:r_i,1:rp] %*% solve(EZZ_BB[1:rp,1:rp])
    Q_i[1:r_i,1:r_i] <- (EZZ[1:r_i,1:r_i] - A_i[1:r_i,1:rp] %*% t(matrix(EZZ_FB[1:r_i,1:rp],r_i,rp))) / TT
    
    A_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- A_i #rp1=0
    Q_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- Q_i;
    V_0_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- Vsmooth[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC),1]
    
  }
  
  rp1 <- sum(sum(r))*ppC #20
  niM <- sum(i_idio[1:nM]) #23
  
  # idiosyncratic
  EZZ <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],2:dim(Vsmooth)[3]], MARGIN = 1:2, FUN = sum))) #E(Z'Z)
  EZZ_BB <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],1:(dim(Vsmooth)[3]-1)], MARGIN = 1:2, FUN = sum)))  #E(Z(-1)'Z_(-1))
  EZZ_FB <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)]))) + diag(diag(apply(VVsmooth[(rp1+1):dim(VVsmooth)[1],(rp1+1):dim(VVsmooth)[2],], MARGIN = 1:2, FUN = sum))) #E(Z'Z_(-1)) 
  
  A_i <- EZZ_FB %*% diag(1/diag(EZZ_BB)) 
  Q_i <- (EZZ - A_i %*% t(EZZ_FB)) / TT 
  
  A_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- A_i[1:niM,1:niM] #更新A[21:43,21:43]
  Q_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- Q_i[1:niM,1:niM]
  V_0_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- diag(diag(Vsmooth[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM),1]))
  
  Z_0 <- Zsmooth[,1] 
  
  nanY<-is.na(y)
  y[nanY] <- 0 
  
  # LOADINGS
  C_new <- C
  
  # Blocks
  bl <- unique(blocks) #去除重复的行。
  n_bl <- dim(bl)[1] #4
  bl_idxM <- NULL
  bl_idxQ <- NULL
  R_con <- NULL
  q_con <- NULL
  
  for(i in 1:n_b){
    bl_idxQ <- cbind(bl_idxQ, matlab::repmat(bl[,i],1,r[i]*ppC)) #bl_idxQ 为 1,1,1（实证）
    bl_idxM <- cbind(bl_idxM, matlab::repmat(bl[,i],1,r[i]), matlab::zeros(n_bl,r[i]*(ppC-1)))
    if(i == 1){
      R_con <- kronecker(R_mat,diag(r[i]))
    }else{
      R_con <- as.matrix(Matrix::bdiag(R_con, kronecker(R_mat,diag(r[i]))))
    }
    q_con <- rbind(q_con, matlab::zeros(r[i]*dim(R_mat)[1],1))
  }
  
  bl_idxM <- bl_idxM == 1  
  bl_idxQ <- bl_idxQ == 1  
  
  #idio
  i_idio_M <- i_idio[1:nM] #23个T
  n_idio_M <- length(which(i_idio_M)) #23
  c_i_idio <- cumsum(i_idio)  
  
  for(i in 1:n_bl){ 
    bl_i <- bl[i,]  
    rs <- sum(r[bl_i == 1]) 
    
    idx_i <- NULL
    for(k in 1:nrow(blocks)){
      idx_i[k] <- sum(blocks[k,] == bl_i) == dim(blocks)[2]
    }
    idx_i <- which(idx_i)  
    
    # MONTHLY--------------
    idx_iM <- idx_i[idx_i < (nM + 1)] 
    n_i <- length(idx_iM)#3
    
    if(n_i==0){
      denom <- NULL
      nom <- NULL
    } else {
      denom <- matlab::zeros(n_i*rs,n_i*rs) #6*6
      nom <- matlab::zeros(n_i,rs) #3*2
      
      
      i_idio_i <- i_idio_M[idx_iM] == 1 #3个T
      i_idio_ii <- c_i_idio[idx_iM]
      i_idio_ii <- i_idio_ii[i_idio_i] 
      
      temp = which(c(bl_idxM[i,],rep(0,nrow(Zsmooth)-ncol(bl_idxM)))==1)
      for(t in 1:TT){ #TT=385
        # t=1
        nanYt <- diag(!nanY[idx_iM,t]) 
        denom <- denom + kronecker(Zsmooth[temp,t+1] %*% t(Zsmooth[temp,t+1]) + Vsmooth[temp,temp,t+1], nanYt)
        nom <- nom + y[idx_iM,t] %*% t(Zsmooth[temp,t+1]) - nanYt[,i_idio_i] %*% (Zsmooth[rp1+i_idio_ii,t+1] %*% t(Zsmooth[temp,t+1]) + Vsmooth[rp1+i_idio_ii,temp,t+1])
        
      }
      
      vec_C <- solve(denom) %*% c(nom)
      C_new[idx_iM,temp] <- matrix(vec_C,n_i,rs)
    }
    
    # QUARTERLY------------------ 
    idx_iQ <- idx_i[idx_i>c(nM)] #idx_i为1,2,6,25。idx_iQ=25
    rps <- rs*ppC #10
    
    R_con_i <- R_con[,bl_idxQ[i,]]
    q_con_i <- q_con #16*1的0
    no_c <- !(rowSums(R_con_i == 0) == ncol(R_con_i))
    R_con_i <- R_con_i[no_c,]
    q_con_i <- q_con_i[no_c,]
    
    ####
    for(j in idx_iQ){ #idx_iQ=25
      denom <- matlab::zeros(rps,rps)
      nom <- matlab::zeros(1,rps)
      idx_jQ <- j-nM
      
      
      if(1 == 1){
        i_idio_jQ <- (rp1+n_idio_M+3*(idx_jQ-1)+1):(rp1+n_idio_M+3*idx_jQ)
        #rp1=20,对于j=25,49:53;对于j=24,44:48。39:41（实证）
        V_0_new[i_idio_jQ,i_idio_jQ] <- Vsmooth[i_idio_jQ,i_idio_jQ,1]
        A_new[i_idio_jQ[1],i_idio_jQ[1]] <- A_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]#rp1=20
        Q_new[i_idio_jQ[1],i_idio_jQ[1]] <- Q_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]
        
        temp = which(c(bl_idxQ[i,],rep(0,nrow(Zsmooth)-ncol(bl_idxQ)))==1)#1,2,3（实证）
        for(t in 1:TT){
          nanYt <- as.vector(!nanY[j,t])*1
          # nn2 <- sum(bl_idxQ[i,])
          denom <- denom + kronecker(Zsmooth[temp,t+1] %*% t(Zsmooth[temp,t+1]) + Vsmooth[temp,temp,t+1],nanYt)
          nom <- nom + y[j,t] %*% t(Zsmooth[temp,t+1])
          nom <- nom - nanYt %*% (matrix(c(1,1,1), nrow = 1) %*% Zsmooth[i_idio_jQ,t+1] %*% t(Zsmooth[temp,t+1]) +
                                    matrix(c(1,1,1), nrow = 1) %*% Vsmooth[i_idio_jQ,temp,t+1])
        }
        C_i <- solve(denom) %*% t(nom)
        C_i_constr <- C_i - solve(denom) %*% t(R_con_i) %*% solve(R_con_i %*% solve(denom) %*% t(R_con_i)) %*% (R_con_i %*% C_i - q_con_i)
        C_new[j,temp] <- C_i_constr
      }
    }
  }
  
  #更新R
  R_new <- matlab::zeros(n,n)
  for(t in 1:TT){
    nanYt <- diag(!nanY[,t])*1 == 1
    R_new <- R_new + (y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) %*% t(y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) 
    + nanYt %*% C_new %*% Vsmooth[,,t+1] %*% t(C_new) %*% nanYt + (diag(n)-nanYt) %*% R %*% (diag(n)-nanYt)
  }
  
  R_new <- R_new/TT
  RR <- diag(R_new) #RR(RR<1e-2) = 1e-2;
  RR[i_idio_M] <- 1e-04
  if(nM<length(RR)){
    RR[(nM+1):length(RR)] <- 1e-04
  }
  R_new <- diag(RR)
  
  if(!is.matrix(Z_0)){
    Z_0<-matrix(Z_0,length(Z_0),1)
  }
  
  # output
  return(list(C_new = C_new, R_new = R_new, A_new = A_new, Q_new = Q_new, Z_0 = Z_0, V_0 = V_0_new, loglik = loglik))
}



# 收敛条件
em_converged <- function(loglik = NULL, previous_loglik = NULL, threshold = NULL, check_increased = NULL){
  
  nargin <- 4 - sum(c(is.null(loglik), is.null(previous_loglik), is.null(threshold), is.null(check_increased)))
  
  if(nargin < 3){threshold <- 1e-4}
  if(nargin < 4){check_increased <- 1}
  
  converged <- 0
  decrease <- 0
  
  if(!is.null(check_increased)){
    if(loglik - previous_loglik < -1e-3){ # allow for a little imprecision
      message(paste("******likelihood decreased from",round(previous_loglik,4),"to",round(loglik,4)))
      decrease <- 1
    }
  }
  
  delta_loglik <- abs(loglik - previous_loglik)
  avg_loglik <- (abs(loglik) + abs(previous_loglik) + 2.2204e-16)/2
  
  if((delta_loglik / avg_loglik) < threshold){converged <- 1}
  
  # output
  list(converged = converged, decrease = decrease)
  
}
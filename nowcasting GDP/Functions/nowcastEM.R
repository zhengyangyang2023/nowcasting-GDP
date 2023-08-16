# 实时预测函数-----

nowcastEM <- function(formula, data, r = NULL, q = NULL, p = NULL, method = 'EM', blocks = NULL, frequency = NULL,max_iter= NULL){
  
  library(matlab)
  # Checking user inputs
  
  # check formula
  if(is.character(formula)){
    formula <- as.formula(formula)
  }
  
  # the number of factors, shocks, and lags
  if(is.null(q) | is.null(r) | is.null(p)){
    warnings('Parameters q, r and p must be specified.')
  }
  
  # the frequencies of the variables
  if(length(frequency) != ncol(data)){
    stop("the length of the frequency vector must be the same as the number of variables in the data object")
  }
  
  if(sum(!frequency%in% c(12,4))!=0){
    stop("The frequencies should be a vector of numerics taking values 4 (quarterly) or 12 (monthly)")
  }
  
  # preparing the data
  k <- model.frame(formula, data, na.action = NULL) #25个指标
  
  y_position <- which(colnames(data) == colnames(k)[1])
  freq_y <- frequency[y_position]
  
  if(freq_y == 4){
    y <- month2qtr(ts(k[,1], start = start(data), frequency = 12))
  }else{
    y <- ts(k[,1], start = start(data), frequency = 12)
  }  
  
  
  # selecting the method
  
  if(method == 'EM'){
    
    # checking validity of inputs
    if(p > 5){
      stop('Parameter p must be less or equal to 5.')
    }
    
    if(is.null(blocks)){
      stop("The block structure determining which variables load into which factors should be specified.")
    }
    
    if(!is.null(q)){
      message("Obs: for this estimation method the number of common shocks is assumed to be equal to the number of factors, i.e. q = r.")
    }
    
    blocks <- as.matrix(blocks) # variables should be rows, columns the blocks
    n_blocks <- dim(blocks)[2]
    
    x <- ts(model.frame(formula, data, na.action = NULL), start = start(data), frequency = frequency(data))
    
    # new frequency
    new_frequency <- c(frequency[y_position], frequency[-y_position])
    
    # reshuffling the blocks
    blocks <- rbind(matrix(blocks[y_position,],ncol = n_blocks),matrix(blocks[-y_position,], ncol = n_blocks))
    
    # determine the number of quarterly series
    nQ <- sum(new_frequency==4)
    
    # preparing X
    
    # 1) all quarterly series should be positioned at the last columns
    # reshuffle vector
    idx <- cumsum(rep(1,dim(x)[2]))
    V_Q <- which(new_frequency==4)
    if(is.integer(V_Q) && length(V_Q) == 0){
      idx_new <- idx
    }else{
      idx_M <- idx[-V_Q]
      idx_new <- c(idx_M,V_Q) #指标索引，月度在前，季度在后
    }
    
    # adapting data base
    x <- x[,idx_new]
    # adapting blocks
    blocks <- as.matrix(blocks[idx_new,]) #将季度频率指标的blocks行放在最后几行
    
    # 2) position of the target variable
    y_pos <- which(colnames(x)==colnames(k)[1])
    Par <- list(r = rep(r,n_blocks), # Number of common factors
                p = p, # Number of lags in autoregressive of factor (same for all factors)
                max_iter = max_iter, # max number of itereations for the EM loop
                i_idio = c(rep(T,dim(x)[2]-nQ), rep(F,nQ)),
                Rconstr = matrix(c(1,1,-1,0,0,-1),2,3), 
                q = matrix(rep(0,2),2,1), 
                nQ = nQ, # Number of quarterly series
                blocks = blocks # Block loadings
    )
    Res <- EM_DFM_SS_block_idioQARMA_restrMQ(x,Par)
    # recovering the factors
    idx_factor <- c(1:r)
    if(dim(blocks)[2]>1){
      idx_factor_aux <- lapply(X = seq(1,dim(blocks)[2]-1), FUN = function(x){(x*r*5 + 1):(x*r*5 + r)})
      for(j in 1:length(idx_factor_aux)){idx_factor <- append(idx_factor, idx_factor_aux[[j]])}
    }
    
    # Factors and estimated explanatory variables
    factors <- list(dynamic_factors = ts(Res$FF[,idx_factor], start = start(x), frequency = 12))#397*4,idx_factor为1,6,11,16
    names(factors$dynamic_factors) <- as.vector(sapply(X = 1:dim(blocks)[2],FUN = function(X){paste0("Block",X,"_factor",1:r)}))
    
    fore_x <- ts(Res$X_sm, start = start(x), frequency = 12)#397*25，后面的缺失都已预测好
    colnames(fore_x) <- colnames(x)
    # Fitted Values
    # finding varepsilon in the FF matrix
    ncol_eps <- nQ*2+ncol(x) 
    idx_eps <- (ncol(Res$FF)-ncol_eps+1):ncol(Res$FF)
    idx_eps_y <- ncol(Res$FF)-ncol_eps+which(Res$C[y_pos,idx_eps]>0)
    
    # finding the autocorelation coefficient of the AR(1) process imposed on varepsilon
    alpha <- solve(t(Res$FF[-nrow(Res$FF),min(idx_eps_y)])%*%Res$FF[-nrow(Res$FF),min(idx_eps_y)])%*%(t(Res$FF[-nrow(Res$FF),min(idx_eps_y)])%*%Res$FF[-1,min(idx_eps_y)])
    
    # y monthly
    if(new_frequency[idx_new[y_pos]]==12){
      # expected error in next period
      error_hat <- c(0,rep(alpha,nrow(Res$FF)-1)*Res$FF[-nrow(Res$FF),min(idx_eps_y)])
      yprev <- ts(Res$FF[,1:(ncol(Res$FF)-ncol_eps)]%*%Res$C[y_pos,1:(ncol(Res$FF)-ncol_eps)]+error_hat, start = start(x), frequency = 12)
      
      
      # denormalize the fitted values
      yprev <- yprev*Res$Wx[y_pos]+Res$Mx[y_pos] #去除标准化 
      yprevmonth <- c()
      yprevmonth1 <- c() 
      # observed values
      y <- x[,y_pos]
    }
    
    # y quarterly
    if(new_frequency[idx_new[y_pos]]==4){
      yprev <- ts(as.vector(Res$FF[,1:(ncol(Res$FF)-ncol_eps)]%*%Res$C[y_pos,1:(ncol(Res$FF)-ncol_eps)]+rowSums(Res$FF[,idx_eps_y])), start = start(x), frequency = 12)
      
      yprev1 <- ts(as.vector(as.matrix(Res$FF[,1:r],ncol=r)%*%Res$C[y_pos,1:r]+Res$FF[,idx_eps_y[1]]), start = start(x), frequency = 12)
      yprevmonth1 <- yprev1*Res$Wx[y_pos]*3+Res$Mx[y_pos] #月度GDP同比
      
      # denormalize the fitted values
      yprev <- yprev*Res$Wx[y_pos]+Res$Mx[y_pos]
      yprevmonth <- yprev[c(length(yprev)-2,length(yprev)-1,length(yprev))]
      # quarterly values
      yprev <- month2qtr(yprev)#nowcasting::month2qtr(yprev,reference_month = "mean")
      
      # observed values
      y <- month2qtr(x[,y_pos])
    }
    
    # Observed and forecast y
    Y <- cbind(y,yprev,yprev)
    Y[is.na(Y[,1]),2] <- NA
    Y[!is.na(Y[,1]),3] <- NA
    colnames(Y) <- c('y','in','out')
    
    # forecast X
    fore_x <- fore_x[,-y_pos]
    colnames(fore_x) <- colnames(x)[-y_pos]
    
    res <- list(yfcst = Y, 
                yfcstmonth = yprevmonth,
                yfcstmonth1 = yprevmonth1,
                factors = factors, 
                xfcst = fore_x,
                Res = Res
    )
    
  }
  
  # output
  return(res)
  
}

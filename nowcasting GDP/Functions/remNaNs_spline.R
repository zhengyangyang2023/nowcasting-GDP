remNaNs_spline <-function(X,options){
  
  #options=optNaN
  #X=x
  
  TT <- dim(X)[1]
  N <- dim(X)[2]
  k <- options$k
  indNaN <- is.na(X)
  
  if(options$method == 1){ # replace all the missing values 
    for (i in 1:N){  
      x = X[,i]
      x[indNaN[,i]] = median(x,na.rm = T);
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
    
  }else if(options$method == 2){ # replace missing values after removing leading and closing zeros
    
    rem1 <- (rowSums(indNaN)>N*0.8) 
    nanLead <- which(rem1) 
    nanLE<-nanLead
    if(length(nanLE)>0){
      X<-X[-nanLE,] 
    }
    
    indNaN=is.na(X)
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx)) 
      t2 = max(which(!isnanx)) 
      
      x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1)) #将第一个非缺失值到最后一个非缺失值的中间缺失值补全
      # xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T) #两头的缺失值使用中位数填充
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      #将x前后各填充k个值，然后使用移动平均。
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]] #原数据的缺失值 使用 移动平均后相应位置的元素 代替。
      X[,i]=x;
    }
    
  }else if(options$method == 3){ # only remove rows with leading and closing zeros
    
    rem1 <- (rowSums(indNaN)==N) # indNaN为原样本缺失情况。N=25
    nanLead <- which(rem1) #每一行全部缺失的行索引
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    if(length(nanLE) != 0){
      X <- X[-nanLE,] #去除掉全部缺失的行构成的样本
    }
    indNaN=is.na(X) #去除掉全部缺失的行后构成的样本 的缺失情况
    
  }else if(options$method == 4){ # remove rows with leading and closing zeros & replace missing values
    
    rem1 <- (rowSums(indNaN)==N)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    X<-X[-nanLE,]
    indNaN=is.na(X)
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
    
  }else if(options$method == 5){
    indNaN=is.na(X)
    
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
  }
  
  return(list(X = X,indNaN=indNaN))
  
}
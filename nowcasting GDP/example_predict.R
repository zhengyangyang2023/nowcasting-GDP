
# 预测---------------
library(openxlsx)
dataindex <- read.xlsx("./Input/model_spec.xlsx")

process_fun <- function(x, method){
  n = length(x)
  if(method == 'ch1'){
    res = c(NA,x[2:n]/x[1:(n-1)]-1)
  }else if(method == 'lin'){
    res = x
  }
  return(res)
}

process_data <- function(date1){
  data0 <- read.xlsx(paste0('./Input/',date1,'.xlsx'),1,rowNames = T,detectDates = T)
  name1 = names(data0)
  data1 <- data0
  block1 = matrix(NA,nrow = ncol(data0),ncol = 4)
  for(i in 1:ncol(data0)){
    temp = dataindex[dataindex['SeriesID'] == name1[i],]
    data1[,i] = process_fun(data0[,i], method = temp[,'Transformation'])
    block1[i,] <- as.numeric(temp[,c("Block1-Global","Block2-Soft","Block3-Real","Block4-Price_Currency")])
  }
  res <- list()
  res[['process_data']] = data1
  res[['process_block']] = block1
  
  return(res)
}

mynowcast <- function(data0,h){ #h为向前预测的期数，h=0预测当前季度，h=1预测当前季度和下一个季度
  mod1 = nrow(data0)%%3
  if(mod1==0){mod1 = 3}
  
  k = ncol(data0)-1 #月度指标个数
  # 在后面填充NA行，最终数据的行数必须为3的倍数
  NA1 = as.data.frame(matrix(NA, ncol=k+1, nrow=3*h+3-mod1))
  names(NA1) <- names(data0)
  data1 <- rbind(data0,NA1)
  data2 <- ts(data1,frequency = 12,start = c(2009,1))
  frequency1 = c(4,rep(12,k))
  
  nowEM_tongbi <- nowcastEM(formula = GDP ~ ., data = data2, r = 1, p = 1, 
                            method = "EM", blocks = block1, frequency = frequency1,max_iter=500)
  pre <- nowEM_tongbi[["yfcst"]]
  time1 = zoo::as.Date(pre)[(nrow(pre)-h):nrow(pre)]
  date1 = zoo::as.yearqtr(time1,format = "%Y-%m-%d")
  nowpre <- data.frame(date = date1, forecast = pre[(nrow(pre)-h):nrow(pre),3])
  return(nowpre)
}

date0 = read.table('./Input/date.txt') #数据截止日期
h0 = read.table('./Input/horizon.txt') #向前预测期数
h = h0[1,1]
for(date1 in date0[,1]){
  res1 = process_data(date1)
  data1 = res1$process_data
  block1 = res1$process_block
  res = mynowcast(data1,h)
  write.xlsx(res,paste0('./Output/',date1,'_horizon_',h,'.xlsx'),rowNames = FALSE)
  message(date1,' has been output')
}



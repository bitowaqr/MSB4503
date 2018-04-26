
# # # # # # # # # # # # # # # # # # # # # # # # # # 
# Using Ridge Regression                          #
# to draw the Training and Test Error curves      # 
# as a function of the amount of regularization   #
# # # # # # # # # # # # # # # # # # # # # # # # # # 


# load packages
  library(caret)
  library(MASS)
  library (glmnet)
  library(doParallel)
  library(cowplot)
  library(glmnet)
  


# parallel
  no_cores <- detectCores() - 1  
  cl <- makeCluster(no_cores, type="FORK")  
  registerDoParallel(cl) 

# load data
  Y = read.table("https://web.stanford.edu/~hastie/Papers/LARS/diabetes.data",header = T)$Y
  data = read.table("https://web.stanford.edu/~hastie/Papers/LARS/data64.txt",header = T)
  data$Y = as.numeric(scale(Y))


#   # plot.list = list()
#   # for(p in 1:6){   # looping over 6 different training/test sets


# split into training and test
  y = names(data) == "Y"
  train.sample.index = sample(x=1:length(data[,1]),
                              size = length(data[,1]) * 2/3) # 2/3 of the data is used as test data
  
  df.train = data[train.sample.index,]
  x.train = df.train[,!y]
  x.train = as.matrix(x.train) # ridge likes matix form...
  
  df.test = data[-train.sample.index,]
  x.test = df.test[,!y]
  x.test = as.matrix(x.test)
  
  y.train = df.train[,y]
  y.test = df.test[,y]

# Fit models:
  fit.ridge <- glmnet(x.train, 
                      y.train, 
                      family="gaussian", 
                      alpha = 0) # = ridge 
# plot all betas
  plot(fit.ridge)
# plot individual beta
  # (probably a unreasonably complicated way to do this...)
  rownames(fit.ridge$beta)
  look.at = c("bmi","age","age.bmi")
  select.var = which( rownames(fit.ridge$beta) %in% look.at)
  var.selected = fit.ridge$beta[select.var,]
  var.selected = data.frame(apply(var.selected,1,function(x) x))
  var.selected$s = rownames(var.selected)
  var.selected = reshape2::melt(data.frame(var.selected))
  var.selected$s = as.numeric(gsub("s","",var.selected$s))
  
  ggplot(var.selected) + 
    geom_line(aes(x=s,y=value,col=variable)) +
    geom_hline(yintercept = 0) + 
    theme_minimal()

# Set lambda grid for CV
  lambda.seq = c(seq(0.0001, 1, length = 100))

# 5-fold Cross validation 
  fit.ridge.cv = train(x = x.train,
                        y = y.train,
                        method = "glmnet", 
                        trControl = trainControl(method="cv",number=5,allowParallel = T),
                        tuneGrid = expand.grid(lambda = lambda.seq, alpha=0)
                        )
  
# retrive train, test and cv error
  RMSE.plot = data.frame(lambda=NA,RMSE.train=NA,RMSE.test = NA)
  for (i in lambda.seq) {
    temp.pred.train = predict(fit.ridge,newx=x.train,s=i)
    temp.rmse.train = sqrt(mean((temp.pred.train - y.train)^2))
    
    temp.pred.test = predict(fit.ridge,newx=x.test,s=i)
    temp.rmse.test = sqrt(mean((temp.pred.test - y.test)^2))
    RMSE.plot = rbind(RMSE.plot,
                      data.frame(lambda=i,
                                 RMSE.train= temp.rmse.train,
                                 RMSE.test = temp.rmse.test))
    
    
  }
  RMSE.plot = RMSE.plot[-1,]
  cv.error = fit.ridge.cv$results$RMSE
  cv.lambda = fit.ridge.cv$results$lambda
  
# plot results
  p = 
    ggplot() +
      geom_line(aes(x=cv.lambda,y=cv.error,col="CV error")) +
      geom_line(data = RMSE.plot, aes(x=lambda,y=RMSE.train,col="Training error")) +
      geom_line(data = RMSE.plot, aes(x=lambda,y=RMSE.test,col="Test error")) +
      scale_x_reverse()

  p

# 
# 
# plot.list[[p]] = ggplot() +
#                     geom_line(aes(x=cv.lambda,y=cv.error,col="CV error")) +
#                     geom_line(data = RMSE.plot, aes(x=lambda,y=RMSE.train,col="Training error")) +
#                     geom_line(data = RMSE.plot, aes(x=lambda,y=RMSE.test,col="Test error")) +
#                     scale_x_reverse()
# }
# plot_grid(plotlist = plot.list,ncol=2)
  
  
  
  
  
  
  # building my own ridge function, using the smaller data set
  easy.dat  = read.table("https://web.stanford.edu/~hastie/Papers/LARS/diabetes.data",header = T)
  easy.dat = data.frame(apply(easy.dat,2,function(x){scale(x,center=T,scale=T)})) # scale and center
  Y = easy.dat$Y # get Y
  easy.dat = easy.dat[,-length(easy.dat)] # remove Y from X matrix
  X = as.matrix(easy.dat, ncol=dim(easy.dat)[2])  #make it a matrix
  
  my.ridge = function(X,y,lambda){
    beta.hat = solve(crossprod(X) + lambda * diag(ncol(X))) %*% crossprod(X,Y)
    y.hat = rowSums(beta.hat * data.frame(X))
    list(beta.hat=beta.hat,
         y.hat=y.hat)
  }
  
  results = my.ridge(X,Y,lambda = 1)
  
  # using the glmnet ridge function:
  beta.hat.glmnet = glmnet(x=X,y=Y,alpha=0,lambda = 1, intercept = F)$beta
  
  # WHY? whats the difference ???
  cbind(my.results=results$beta.hat,
        glmnet=beta.hat.glmnet)
  

  


####### GRADIENT DECENT
# set a random start
# use gradient decent to improve slope and intercept
# plot the the estimate every 50 iterations
# plot MSE over all iterations

gradientDesc <- function(x=mtcars$disp, 
                         y=mtcars$mpg, 
                         learn_rate = 0.00293, 
                         conv_threshold =  0.001, 
                         max_iter=2500000) {
  MS.save = numeric()
  n = length(x)
  plot(x, y, col = "blue", pch = 20,xlim=c(-5,5),ylim=c(-5,5))
  plot.iter = 0
  m <- runif(1, -5, 5)
  c <- runif(1, -5, 5)
  abline(c, m) 
  yhat <- m * x + c
  MSE <- sum((y - yhat) ^ 2) / n
  converged = F
  iterations = 0
  while(converged == F) {
    ## Implement the gradient descent algorithm
    plot.iter = plot.iter + 1
    m_new <- m - learn_rate * ((1 / n) * (sum((yhat - y) *x)))
    c_new <- c - learn_rate * ((1 / n) * (sum(yhat - y)))
    m <- m_new
    c <- c_new
    if(plot.iter==50){
      plot.iter = 0
      abline(c, m) 
    }
    yhat <- m * x + c
    MSE_new <- sum((y - yhat) ^ 2) / n
    MS.save[iterations] = MSE_new
    if(MSE - MSE_new <= conv_threshold) {
      abline(c, m) 
      converged = T
      print(paste("Optimal intercept:", c, "Optimal slope:", m))
      return(MS.save)
    }
    iterations = iterations + 1
    if(iterations > max_iter) { 
      abline(c, m) 
      converged = T
      print(paste("Optimal intercept:", c, "Optimal slope:", m))
      return(MS.save)
    }
  }
}


# Run the function 
learn_rate = 0.00293
conv_threshold =  0.001
max_iter=10000

par(mfrow=c(3,2))
for(i in 1:3){
  plot(
    gradientDesc(scale(mtcars$disp), 
                 scale(mtcars$mpg), 
                 learn_rate,
                 conv_threshold,
                 max_iter),
    type="l",
    log="y",
    ylab="log(MSE)",
    xlab="iteration",
    main="MSE per iteration")
}


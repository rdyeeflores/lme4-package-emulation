
emu.lmer <- function(outcome, predictor, group, data){

  outcome <- substitute(outcome)
  predictor <- substitute(predictor)
  group <- substitute(group)  
  data <- data[complete.cases(data), ]; data <- data[order(data[,as.character(group)],decreasing=FALSE),]
  
  J <- length(unique(data[,as.character(group)]))
  I <- as.numeric(table(data[,as.character(group)])[1])

  ## Packages
  library(mvtnorm) 
  library(numDeriv)
  
  ## Data partitioning   
  y <- data[,as.character(outcome)] 
  X <- cbind(rep(1, length(y)), data[,as.character(predictor)])

  model_negloglikelihood <- function(param){
    
    ## Getting Z, covariance matrix of random effects, & G
    Z <- matrix(0, nrow=I*J, ncol=2*J)
    for (i in 1:J) {
      Z[(I*i-I+1):(I*i),(i*2-1):(i*2)] <- X[(I*i-I+1):(I*i),]}
    CVM.rand <- cbind(c(param[4], param[6]), c(param[6], param[5]))
    G <- matrix(0, nrow=2*J, ncol=2*J)
    for (i in 1:J) {
      G[(i*2-1):(i*2),(i*2-1):(i*2)] <- CVM.rand}
    
    ## Getting negative log-likelihood using BETA, Z, G, and r.vari
    ZGZt <- Z %*% G %*% t(Z)
    XB <- X %*% param[1:2]
    irv <- diag(I*J) * param[3]
    LL <- dmvnorm(t(y), mean = XB, sigma = ZGZt + irv, log=TRUE)
    (negLL <- -LL)
    
  }
  
  ## Initial guesses for nlminb() 
  BETA <- solve(crossprod(X)) %*% t(X) %*% y    
  r <- y - X %*% BETA                          
  r.vari <- sum(r^2) / (length(y) - 2)   
  var_inter <- 0
  var_slope <- 0
  covari <- 0
  
  ## Running nlminb() 
  OPTIM <- nlminb(start = c(BETA, r.vari, var_inter, var_slope, covari), 
                  objective = model_negloglikelihood, 
                  control = list(eval.max=1e3, iter.max=1e3)) ## both at about 5x
  
  #### TABLE
  
  ## Random effects: variances, SDs, and correlation
  rev <- OPTIM$par[4:5]
  resd <- sqrt(OPTIM$par[4:5])
  rec <- OPTIM$par[6] / sqrt(OPTIM$par[4]*OPTIM$par[5])
  ## Residual variance and SD
  rv <- OPTIM$par[3]
  rsd <- sqrt(OPTIM$par[3])
  ## Fixed effects: estimates, SEs, and t-values
  fee <- OPTIM$par[1:2]
  H <- hessian(func=model_negloglikelihood, x=OPTIM$par)
  fese <- sqrt(diag(solve(H)))[1:2]
  fet <- OPTIM$par[1:2] / sqrt(diag(solve(H)))[1:2]
  ## Correlation of fixed effects
  fec <- solve(H)[2,1] / sqrt(solve(H)[1,1]*solve(H)[2,2])
  
  re <- rbind(c(rev[1], resd[1], NA), c(rev[2], resd[2], rec), c(rv, rsd, NA))
  rownames(re) <- c("(Intercept)", "Predictor", "Residual")
  colnames(re) <- c("Variance", "SD", "Correlation")
  
  fe <- rbind(c(fee[1], fese[1], fet[1], NA), c(fee[2], fese[2], fet[2], fec))
  rownames(fe) <- c("(Intercept)", "Predictor")
  colnames(fe) <- c("Estimate", "SE", "t-value", "Correlation")
  
  list("Random effects", re, "Fixed effects", fe)

}


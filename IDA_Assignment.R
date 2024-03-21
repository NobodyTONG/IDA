##### Huantong Hou (s2481591)

##### Problem 2
## load data and library
load("/Users/houhuantong/Edinburgh/S2/IDA/Assignment/dataex2.Rdata")
library(mice)

## set necessary parameters
beta0=1
beta1=3
nsim = 100
n = 100
seed = 1
M=20

## Stochastic Regression Model
# Extract X and Y from each subset
X <- Y <- matrix(0, nrow=n, ncol=nsim)
for(l in 1:nsim){
  X[, l] <- dataex2[, 1, l]
  Y[, l] <- dataex2[, 2, l]
}
estimates_SRI <- matrix(0, nrow = 4, ncol = nsim)
for(l in 1:nsim){
  # extract columns without NA
  nonNA <- which(!is.na(dataex2[, 2, l]))
  # fit data excluding NA
  fit <- lm(Y[nonNA, l] ~ X[nonNA, l])
  sigmaest <- sigma(fit)
  # predict values using fitting model
  pred <- fit$coefficients[1] + dataex2[, 1, l]*fit$coefficients[2] + rnorm(n, 0, sigmaest)
  dataex2_completed <- Y[, l]
  Y_NA <- which(is.na(dataex2[, 2, l]))
  dataex2_completed[Y_NA] <- pred
  estimates_SRI[1, l] <- fit$coefficients[2]
  estimates_SRI[2, l] <- var(dataex2_completed)
  estimates_SRI[3, l] <- estimates_SRI[1, l] - qt(0.975, n-1)*sqrt(estimates_SRI[2, l]/n)
  estimates_SRI[4, l] <- estimates_SRI[1, l] + qt(0.975, n-1)*sqrt(estimates_SRI[2, l]/n)
}
# calculate empirical coverage probability of SRI
coverage_SRI <- sum((estimates_SRI[3,] <= beta1) & (estimates_SRI[4,] >= beta1))/nsim

## Bootstrap
estimates_boot <- matrix(0, nrow = 2, ncol = nsim)
# for each data set
for (i in 1:nsim) {
  df <- data.frame(dataex2[, , i])
  # use Bootstrap to impute missing values
  imp <- mice(df, m = M, meth=c("", "norm.boot"), seed = 1, printFlag = FALSE)
  fits <- with(imp, lm(Y~X))
  ests <- pool(fits)
  # extract lower bound
  estimates_boot[1, i] <- summary(ests, conf.int=TRUE)[2, 7]
  # extract upper bound
  estimates_boot[2, i] <- summary(ests, conf.int=TRUE)[2, 8]
}
# calculate empirical coverage probability of Bootstrap method
coverage_boot <- sum((estimates_boot[1,] <= beta1) & (estimates_boot[2,] >= beta1))/nsim

## Output the results as a table
estimates <- data.frame("SRI" = coverage_SRI,
                        "Bootstrap" = coverage_boot)
rownames(estimates) <- c("Empirical Coverage Probability")
knitr::kable(estimates, escape = FALSE, caption = "Coverage differences")



##### Problem 3 - solutions_workshop_3
load("/Users/houhuantong/Edinburgh/S2/IDA/Assignment/dataex3.Rdata")
require(maxLik)
## Define log likelihood function
log_like <- function(mu, data){
  # extract X and R
  X <- data[, 1]
  R <- data[, 2]
  sigma <- 1.5
  # define probability density function value
  phi <- dnorm(X, mean = mu, sd = sigma)
  # define cumulative density function value
  PHI <- pnorm(X, mean = mu, sd = sigma)
  # get the sum
  sum(R*log(phi) + (1-R)*log(PHI))
}
## obtain the maximum likelihood estimate of miu
mle <- maxLik(logLik = log_like, data = dataex3, start=0)
summary(mle)



##### Problem 5
load("/Users/houhuantong/Edinburgh/S2/IDA/Assignment/dataex5.Rdata")
missing_indicator <- ifelse(is.na(dataex5[,2]), 0, 1)
# EM Algorithm
EM_algorithm <- function(Y_obs, X_obs, start, max_iter = 10000, tol = 1e-7) {
  n <- length(Y_obs)
  beta <- start
  iter <- 1
  converged <- FALSE
  # Loop until the function converges
  while (iter <= max_iter && !converged) {
    # E step
    log_likelihood <- function(beta, Y_obs, X_obs) {
      eta <- beta[1] + X_obs * beta[2]
      p <- exp(eta) / (1 + exp(eta))
      log_likelihood <- sum(Y_obs * log(p) + (1 - Y_obs) * log(1 - p))
      return(-log_likelihood)  
    }
    # M step - maximize the log likelihood function
    result <- optim(par = beta, fn = log_likelihood, Y_obs = Y_obs, X_obs = X_obs)
    new_beta <- result$par
    # justify whether the function is converge
    if (max(abs(beta - new_beta)) < tol) {
      converged <- TRUE
    }
    beta <- new_beta
    iter <- iter + 1
  }
  return(beta)
}
# Calculate the estimate values
Y_obs <- dataex5[which(missing_indicator==1), "Y"]
X_obs <- dataex5[which(missing_indicator==1), "X"]
beta_hat <- EM_algorithm(Y_obs, X_obs, start = c(0, 0))
# Output the results
beta_hat

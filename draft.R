
library(survival)
library(parallel)
library(randomForestSRC)
library(polspline)
library(flexsurv)
library(actuar)

# function simulating data
simulate = function(scenario = 1, n = 5000, t0 = 0.6){
  n.original=n
  n = ceiling(n*1.01)
  if (scenario == 1){
    eps_1 = runif(n,-1,1)
    eps_2 = runif(n)
    eps_3 = runif(n)
    L = 0.5*eps_1 + rnorm(n,0, sd = 1)  
    U = eps_1^2 + rnorm(n,0, sd = 0.3)  
    Z = cbind(L,U)
    t = -log(eps_2)/exp(as.vector(Z%*%c(-1, 2)))
    c = -log(eps_3)/exp(as.vector(Z%*%c(1, 1.5)))
  }
  
  X = pmin(t,c)
  Delta = as.numeric(t < c)
  
  # create two extra columns Tstar = min(T, t0), and Gamma = I(Tstar < C)
  Tstar = pmin(X, t0) # will only use entries that X[i] is event time, X[i] is censor time and the censor time is greater than t0
  Gamma = ifelse(Delta == 1 | (Delta == 0 & t0 < X), 1, 0) 
  

  duplicated_index = unique(c(which(duplicated(X)), which(duplicated(t)), which(duplicated(c))))
  if (length(duplicated_index) != 0) {
    cat(paste(length(duplicated_index) ,'duplicates detected.\n'))
    return(cbind(t, c, X, Delta, L, U)[-duplicated_index,][1:n.original,])
  } else {
    return(cbind(t, c, X, Delta, L, U, Tstar, Gamma)[1:n.original,])
  }
  
}

# create practitioner observed censored data with maximum follow up time tau 
data = simulate(1, 5000)[,-(1:2)]
tau = max(data[,1]) + 1  # set maximum follow-up time
#tau = 1 
data[,2] = data[,2]*(data[,1] < tau)
data[,1] = ifelse(data[,1] < tau, data[,1], tau)



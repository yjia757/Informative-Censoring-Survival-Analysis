# this file is for debugging Cox, Spline, and RSF implemented in sim_full_debug.R
# which were originally used in Jiyu's code to generate matrices S_t, S_c with 
# cross fitting.

# RSF require memory that my PC is not support. So focus on Cox and Spline. 



library(survival)
library(parallel)
library(randomForestSRC)
library(polspline)
library(flexsurv)
library(actuar)

# function simulating data
simulate = function(scenario = 1, n = 1000){
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
  
  duplicated_index = unique(c(which(duplicated(X)), which(duplicated(t)), which(duplicated(c))))
  if (length(duplicated_index) != 0) {
    cat(paste(length(duplicated_index) ,'duplicates detected.\n'))
    return(cbind(t, c, X, Delta, L, U)[-duplicated_index,][1:n.original,])
  } else 
    return(cbind(t, c, X, Delta, L, U)[1:n.original,])
}

# create practitioner observed censored data with maximum follow up time tau 
set.seed(2008)
data = simulate(1, 10000)[,-(1:2)] # Use larger sample size to be on the safe side
tau = max(data[,1]) + 1  # set maximum follow-up time
#tau = 1 
data[,2] = data[,2]*(data[,1] < tau)
data[,1] = ifelse(data[,1] < tau, data[,1], tau)


# function estimate two nuisances S_t and S_c using Cox, Spline, or RSF
nuis_estm = function(data, T_LU = 'Cox', C_LU ='Cox', tau = 1, k = 2){
  # inputs:
  # data: nx4 matrix with columns (X, Delta, L, U) 
  # X - nx1, the censored event time for each patient
  # Delta - nx1, the event indicator for each patient, 1 = Died, 0 = Censored.
  # L - nx1, observed covariate
  # U - nx1, latent covariate
  # T_LU - working model for the T|L,U, available ones include ('Cox', 'Spline', 'RSF')
  # C_LU - working model for the C|L,U, available ones include ('Cox', 'Spline', 'RSF')
  # k - number of folds for cross fitting
  
  data = as.data.frame(data)
  colnames(data)[1:4] = c('X','Delta','L', 'U')
  X.full = data[,1]
  X.sort.full = sort(unique(X.full))
  folds = cut(1:nrow(data), breaks = k, labels = F)
  S_t.temp = S_c.temp = matrix(0, nrow = nrow(data), ncol = length(unique(X.full)))
  
  samp.split_S_t = function(dfs, dfn, X, Delta, L, U, X.sort, X.full, X.sort.full){
    S_t = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
    if (T_LU == 'Cox'){
      model = coxph(Surv(X, Delta) ~ ., timefix = F, data = dfn)
      Lambda = exp(as.matrix(dfs[,-(1:2)])%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
      S_t[,X.sort.full %in% X.sort] = exp(-Lambda)
    }
    if (T_LU == 'Spline'){
      model = hare(X, Delta, cbind(L,U))
      S_t = 1 - sapply(X.sort.full, function(x) phare(q = x, cov = dfs[,-(1:2)], fit = model))
      S_t[is.nan(S_t)] = NA
    }
    if (T_LU == 'RSF'){
      model = rfsrc(Surv(X, Delta) ~ ., data = dfn, ntree = 1000, ntime = 0, mtry = 2, splitrule = 'bs.gradient')
      S_t[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs)$survival
    }
    S_t[,1] = ifelse(S_t[,1] == 0 | is.na(S_t[,1]), 1, S_t[,1]) 
    S_t[S_t==0] = NA
    S_t = t(as.matrix(tidyr::fill(data.frame(t(S_t)), names(data.frame(t(S_t))))))
    return(S_t)
  }
  
  samp.split_S_c = function(dfs, dfn.c, X, Delta_c, L, U, X.sort, X.full, X.sort.full){
    S_c = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
    if (C_LU == 'Cox'){
      model = coxph(Surv(X, Delta_c) ~ ., timefix = F, data = dfn.c)
      Lambda_c = exp(as.matrix(dfs[,-(1:2)])%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
      S_c[,X.sort.full %in% X.sort] = exp(-Lambda_c)
    }
    if (C_LU == 'Spline'){
      model = hare(X, Delta_c, cbind(L,U))
      S_c = 1 - sapply(X.sort.full, function(x) phare(q = x, cov = dfs[,-(1:2)], fit = model))
      S_c[is.nan(S_c)] = NA
    }
    if (C_LU == 'RSF'){
      model = rfsrc(Surv(X, Delta_c) ~ ., data = dfn.c, ntree = 1000, ntime = 0, mtry = 2, splitrule = 'bs.gradient')
      S_c[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs)$survival
    }
    S_c[,1] = ifelse(S_c[,1] == 0 | is.na(S_c[,1]), 1, S_c[,1])
    S_c[S_c==0] = NA
    S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
    return(S_c)
  }
  
  # cross fitting
  for (fold in 1:k){
    dfs = dfs.c = data[folds==fold,]
    dfn = dfn.c = data[folds!=fold,]
    
    X = dfn[,1]
    X.sort = sort(unique(X))
    Delta = dfn[,2]
    Delta_c = (1 - Delta)*(X < tau)
    dfn.c$Delta = Delta_c
    colnames(dfn.c)[2] = 'Delta_c'
    
    X2 = dfs[,2]  # this can be omitted, does not cause difference 
    Delta2 = dfs[,2]  # this can be omitted, does not cause difference
    Delta_c2 = (1 - Delta2)*(X2 < tau)  # this can be omitted, does not cause difference
    dfs.c$Delta = Delta_c2  # this can be omitted, does not cause difference
    colnames(dfs.c)[2] = 'Delta_c'
    
    L = dfn[,3]
    U = dfn[,4]
    
    S_t.temp[folds==fold,] = samp.split_S_t(dfs, dfn, X, Delta, L, U, X.sort, X.full, X.sort.full)
    S_c.temp[folds==fold,] = samp.split_S_c(dfs.c, dfn.c, X, Delta_c, L, U, X.sort, X.full, X.sort.full)
  }
  
  S_t = as.matrix(S_t.temp)
  S_c = as.matrix(S_c.temp)
  
  return(list(S_t, S_c))
}

# estimate of two nuisances 
#nuis = nuis_estm(data, T_LU = 'Spline', C_LU ='Spline', tau, k = 2)
nuis = nuis_estm(data, T_LU = 'Cox', C_LU ='Cox', tau, k = 2) # 5-fold performs better than 2-fold
S_t = nuis[[1]]
S_c = nuis[[2]]

# true nuisances 
# true S_t
tn.X = data[,1]
tn.X.sort = sort(unique(tn.X))
tn.beta.true = c(-1, 2)
tn.Lambda.true = exp(as.matrix(data[,3:4])%*%tn.beta.true) %*% t(tn.X.sort)
tn.surv.true = exp(-tn.Lambda.true)
# true S_c
tn.c.beta.true = c(1, 1.5)
tn.c.Lambda.true = exp(as.matrix(data[,3:4])%*%tn.c.beta.true) %*% t(tn.X.sort)
tn.c.surv.true = exp(-tn.c.Lambda.true)

# compare estimated and truth on observation data[10,]
# S_t
plot(tn.X.sort, tn.surv.true[3,], col = 1, type = "l",lty = 1,
     main = "survival curves",
     xlab = "time", ylab = "Survival probability")
lines(tn.X.sort, S_t[3,], col = 2, type = "l", lty = 2, lwd = 2.5)
legend("topright", 
       legend = c("True Survival of T", "Estimated Cox"),
       col = 1:2, lty = 1:2, lwd = c(1, 2.5), bty = "n")
# S_c
plot(tn.X.sort, tn.c.surv.true[3,], col = 1, type = "l",lty = 1,
     main = "survival curves",
     xlab = "time", ylab = "Survival probability")
lines(tn.X.sort, S_c[3,], col = 2, type = "l", lty = 2, lwd = 2.5)
legend("topright", 
       legend = c("True Survival of C", "Estimated Cox"),
       col = 1:2, lty = 1:2, lwd = c(1, 2.5), bty = "n")








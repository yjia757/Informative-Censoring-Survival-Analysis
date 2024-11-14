# main function for generating data satisfying the latent CAR assumption
# estimate S^2 with cox, spline, or random survival forests for fitting the nuisance parameters 
# conduct sensitivity analysis. 

# focus on estimating theta with both L, U. Once get this right, then move on to only L.

# repeat 50 times with different dataset and average results over for final report 

library(survival)
library(parallel)
library(randomForestSRC)
library(polspline)
library(flexsurv)
library(actuar)

# transformation of the survival time in defining the parameter of interest
nu <- function(t,t0=0.6){
  # # identity function
  # result = t
  
  # indicator function
  result = as.numeric(t<=t0)
  return(result)
}

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
  Delta = as.numeric(t <= c)
  
  # create two extra columns Tstar = min(T, t0), and Gamma = I(Tstar < C)
  Tstar = pmin(X, t0) # will only use entries that X[i] is event time, X[i] is censor time and the censor time is greater than t0
  Gamma = ifelse(Delta == 1 | (Delta == 0 & t0 <= X), 1, 0)  # when Tstar is incorrect, the corresponding Gamma entry is 0
  
  duplicated_index = unique(c(which(duplicated(X)), which(duplicated(t)), which(duplicated(c))))
  if (length(duplicated_index) != 0) {
    cat(paste(length(duplicated_index) ,'duplicates detected.\n'))
    return(cbind(t, c, X, Delta, L, U)[-duplicated_index,][1:n.original,])
  } else {
    return(cbind(t, c, X, Delta, L, U, Tstar, Gamma)[1:n.original,])
  }
}

# find the true parameter of interest using law of large number and simulated event time T
approx_truth = function(n_sim = 500, senario=1, n_obv=20000, t0=0.6){
  store <- numeric(n_sim)
  for (i in 1:n_sim) {
    nu_t <- sapply(simulate(senario, n_obv, t0)[,1], nu, t0)
    store[i] <- mean(nu_t)
  }
  hist_obj <- hist(store, probability = TRUE, 
                   breaks = 50,
                   main = "Histogram", 
                   xlab = "Sample Mean of nu(T)",
                   col = "lightblue", 
                   border = "black")
  abline(v = mean(store), col = "red", lwd = 2)
  text(mean(store), max(hist_obj$density), 
       labels = paste("Mean =", round(mean(store), 4)), 
       pos = 4, col = "red")
  return(mean(store))
}
# approx_truth(n_sim = 100, senario=1, n_obv=20000, t0=1.3)  # 0.638 with t0=0.6 // 0.792 with t0=1.3 //0.908 with t0=3


# function estimate two nuisances S_t and S_c using Cox, Spline, or RSF with full covariates (L,U)
full.nuis_est = function(data, T_LU = 'Cox', C_LU ='Cox', tau = 1, k = 2){
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
  data = data[1:4]
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
    dfs = dfs.c = data[folds==fold,] # Tstar and Gamma are not used in nuisance estimation
    dfn = dfn.c = data[folds!=fold,] # Tstar and Gamma are not used in nuisance estimation
    
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


# function estimate two nuisances S_t and S_c using Cox, Spline, or RSF with only observed covariates L
obv.nuis_est = function(data, T_L = 'Cox', C_L ='Cox', tau = 1, k = 2){
  # inputs:
  # data: nx4 matrix with columns (X, Delta, L, U) 
  # X - nx1, the censored event time for each patient
  # Delta - nx1, the event indicator for each patient, 1 = Died, 0 = Censored.
  # L - nx1, observed covariate
  # T_L - working model for the T|L,U, available ones include ('Cox', 'Spline', 'RSF')
  # C_L - working model for the C|L,U, available ones include ('Cox', 'Spline', 'RSF')
  # k - number of folds for cross fitting
  
  data = as.data.frame(data)
  data = data[1:3]
  colnames(data)[1:3] = c('X','Delta','L')
  X.full = data[,1]
  X.sort.full = sort(unique(X.full))
  folds = cut(1:nrow(data), breaks = k, labels = F)
  S_t.temp = S_c.temp = matrix(0, nrow = nrow(data), ncol = length(unique(X.full)))
  
  samp.split_S_t = function(dfs, dfn, X, Delta, L, X.sort, X.full, X.sort.full){
    S_t = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
    if (T_L == 'Cox'){
      model = coxph(Surv(X, Delta) ~ ., timefix = F, data = dfn)
      Lambda = exp(as.matrix(dfs[,-(1:2)])%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
      S_t[,X.sort.full %in% X.sort] = exp(-Lambda)
    }
    if (T_L == 'Spline'){
      model = hare(X, Delta, L)
      S_t = 1 - sapply(X.sort.full, function(x) phare(q = x, cov = dfs[,-(1:2)], fit = model))
      S_t[is.nan(S_t)] = NA
    }
    if (T_L == 'RSF'){
      model = rfsrc(Surv(X, Delta) ~ ., data = dfn, ntree = 1000, ntime = 0, mtry = 2, splitrule = 'bs.gradient')
      S_t[, which(X.sort.full %in% model$time.interest)] = predict(model, newdata = dfs)$survival
    }
    S_t[,1] = ifelse(S_t[,1] == 0 | is.na(S_t[,1]), 1, S_t[,1]) 
    S_t[S_t==0] = NA
    S_t = t(as.matrix(tidyr::fill(data.frame(t(S_t)), names(data.frame(t(S_t))))))
    return(S_t)
  }
  
  samp.split_S_c = function(dfs, dfn.c, X, Delta_c, L, X.sort, X.full, X.sort.full){
    S_c = data.frame(matrix(0, nrow = nrow(dfs), ncol = length(unique(X.full))))
    if (C_L == 'Cox'){
      model = coxph(Surv(X, Delta_c) ~ ., timefix = F, data = dfn.c)
      Lambda_c = exp(as.matrix(dfs[,-(1:2)])%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
      S_c[,X.sort.full %in% X.sort] = exp(-Lambda_c)
    }
    if (C_L == 'Spline'){
      model = hare(X, Delta_c, L)
      S_c = 1 - sapply(X.sort.full, function(x) phare(q = x, cov = dfs[,-(1:2)], fit = model))
      S_c[is.nan(S_c)] = NA
    }
    if (C_L == 'RSF'){
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
    dfs = dfs.c = data[folds==fold,] # Tstar and Gamma are not used in nuisance estimation
    dfn = dfn.c = data[folds!=fold,] # Tstar and Gamma are not used in nuisance estimation
    
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

    S_t.temp[folds==fold,] = samp.split_S_t(dfs, dfn, X, Delta, L, X.sort, X.full, X.sort.full)
    S_c.temp[folds==fold,] = samp.split_S_c(dfs.c, dfn.c, X, Delta_c, L, X.sort, X.full, X.sort.full)
  }
  
  S_t = as.matrix(S_t.temp)
  S_c = as.matrix(S_c.temp)
  
  return(list(S_t, S_c))
}

# function estimate the parameter of interest theta using estimated nuisance  
para_estm = function(data, S_t.matrix, S_c.matrix, min.S.c = c(1e-4, 1e-3, 1e-2), t0=0.6){
  
  S_t = S_t.matrix
  S_c = S_c.matrix
  
  # define all variables 
  data = as.data.frame(data)
  n = nrow(data)
  X = data[,1]
  X.sort = sort(unique(X))
  n_res = length(unique(X))
  event_rank = pmin(n_res,rank(X, ties.method = 'first'))
  Delta = data[,2]
  Delta_c = (1 - Delta)*(X < tau)
  Y = t(sapply(event_rank, function(x) c(rep(1,x), rep(0,n_res-x))))  # n x n_res matrix  (sample x time)
  Tstar = data[,5]
  Gamma = data[,6]
  col_idx = sapply(Tstar, function(x) {
    return(which.min(abs(x-X.sort)))
  })
  idx = which.min(abs(t0-X.sort))
  nut = matrix(rep(nu(X.sort, t0),n), nrow = n, byrow = TRUE)  # n x n_res  (sample x time)
  
  coll.est_DR = c()
  calc_min.S = function(x) tanh(x)/2+0.5
  for (z in 1:length(min.S.c)) {
    # weight trimming 
    min.S.base = atanh((min.S.c[z] - 0.5)*2) + 0.25
    min.S = calc_min.S(min.S.base)
    if (min(S_t) < min.S) S_t = 1 - (1 - S_t)*(1 - min.S)/(1-min(S_t))
    if (min(S_c) < min.S) S_c = 1 - (1 - S_c)*(1 - min.S)/(1-min(S_c))
    
    # more preparation quantities 
    Lambda_c = -log(S_c)  # n x n_res matrix  (sample x time)
    lambda_c = t(apply(Lambda_c, 1, function(x) c(x[1], diff(x))))  # n x n_res matrix  (sample x time)
    dNc = t(sapply(event_rank, function(x) tabulate(x, nbins = n_res)))*Delta_c  # n x n_res matrix  (sample x time)
    dMc = dNc - Y * lambda_c  # n x n_res matrix  (sample x time)
    dF_t = t(apply(1-S_t, 1, function(x) c(x[1], diff(x))))  # n x n_res  (sample x time)
    
    # more preparation quantities 
    ScTstar = S_c[cbind(1:n, col_idx)]  # n x 1 vector
    # option 1
    S_t.t0 = S_t[,idx]
    J = nut * (S_t - S_t.t0)
    # option 2
    #nut = matrix(rep(nu(X.sort, t0),n), nrow = n, byrow = TRUE)  # n x n_res  (sample x time)
    #value = nut * dF_t
    #J = t(apply(value, 1, function(x) rev(cumsum(rev(x)))))
    
    # compute the denominator of the estimator
    den1 = Gamma/ScTstar  # n x 1 vector 
    den2 = rowSums(dMc/S_c)  # n x 1 value 
    
    Den1 = sum(den1)  # 1 x 1 value 
    Den2 = sum(den2)  # 1 x 1 value
    
    # compute the numerator of the estimator
    num1 = Gamma*nu(X, t0)*Delta/ScTstar  # n x 1 vector
    num2 = rowSums((J/S_t)*(dMc/S_c))  # n x 1 vector
    
    Num1 = sum(num1)
    Num2 = sum(num2)
    
    # estimator based on AIPCW 
    coll.est_DR[z] = (Num1+Num2)/(Den1+Den2)
  }
  
  return(coll.est_DR)
}


# do 100 times and average the results for final report 
covas = c('L,U', 'L')
alg = c('Cox', 'Spline')
chos.t0 = c(0.6)
truth = numeric(length = length(chos.t0))
for (j in 1:length(chos.t0)) {
  truth[j] = approx_truth(n_sim = 100, senario=1, n_obv=20000, t0=chos.t0[j])
}
n_reps = 50
#all.results = vector("list", length(alg))
#for (i in 1:length(alg)) {
#  all.results[[i]] = matrix(0, nrow = length(chos.t0), ncol = n_reps)
#}
summary.table = as.data.frame(matrix(0, nrow = length(covas) * length(alg) * length(chos.t0), ncol = 7))
colnames(summary.table) = c('Covariates', 'Estimator', 't0', 'Truth', 'Esimate', 'Bias', 'SD')

for (a in 1:length(covas)) {

  for (i in 1:length(alg)) { 
    
    for (j in 1:length(chos.t0)) {
      
      results = numeric(n_reps)
      
      for (k in 1:n_reps) {
        set.seed(2000+a+i+j+k)
        # create practitioner observed censored data with maximum follow up time tau 
        data = simulate(scenario = 1, n = 1000, t0 = chos.t0[j])[,-(1:2)]
        tau = max(data[,1]) + 1  # set maximum follow-up time
        #tau = 1 
        data[,2] = data[,2]*(data[,1] < tau)
        data[,1] = ifelse(data[,1] < tau, data[,1], tau)
        
        # estimate of two nuisances 
        if (covas[a] == 'L,U') {
          nuis = full.nuis_est(data, T_LU = alg[i], C_LU = alg[i], tau, k = 2)
        } else if (covas[a] == 'L') {
          nuis = obv.nuis_est(data, T_L = alg[i], C_L = alg[i], tau, k = 2)
        }
        
        S_t = nuis[[1]]
        S_c = nuis[[2]]
        
        # estimate the parameter of interest using the nuisance functions
        coll.est_DR = para_estm(data, S_t.matrix = S_t, S_c.matrix = S_c, 
                                min.S.c = c(1e-3), t0 = chos.t0[j])
        
        est_DR = coll.est_DR[which.min(abs(coll.est_DR-truth[j]))]
        
        results[k] <- est_DR
        
        print(paste("Finished the", k, "th loop of t0 =", chos.t0[j], "and", alg[i], 'with covariates', covas[a]))
      }
      #all.results[[i]][j,] = results
      
      # code for organizing the results
      row_index = length(chos.t0) * length(alg) * (a - 1) + length(chos.t0) * (i - 1) + j
      summary.table[row_index,1] = covas[a]
      summary.table[row_index,2] = alg[i]
      summary.table[row_index,3] = chos.t0[j]
      summary.table[row_index,4] = truth[j]
      summary.table[row_index,5] = mean(results)
      summary.table[row_index,6] = mean(results) - truth[j]
      summary.table[row_index,7] = sd(results)
      
    }
  }
}














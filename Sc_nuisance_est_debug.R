# this file is for debugging Cox, Spline, and RSF implemented in sim_full_debug.R
# which were originally used in Jiyu's code to generate matrices S_c. 

# RSF require memory that my PC is not support. So focus on Cox and Spline. 


library(survival)
library(parallel)
library(randomForestSRC)
library(polspline)
library(flexsurv)
library(actuar)

file_path <- "/Users/yiranjia/Desktop/in_desktop/stats499_lily/Writing/SA/SA_draft1.6/figures"

# function simulating data
simulate = function(scenario = 1, n = 5000){
  n.original=n
  n = ceiling(n*1.01)
  if (scenario == 1){
    eps_1 = runif(n,-1,1)
    eps_2 = runif(n)
    eps_3 = runif(n)
    L = 0.5*eps_1 + rnorm(n,0, sd = 1)  
    U = eps_1^2 + rnorm(n,0, sd = 0.3)  
    Z = cbind(L,U)
    c = -log(eps_2)/exp(as.vector(Z%*%c(-1, 2)))
    t = -log(eps_3)/exp(as.vector(Z%*%c(1, 1.5)))
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


# generate data for debugging purpose, check nuisance estimation
set.seed(2008)
data_debug = simulate(1, 5000)[,-(1:2)] # Use larger sample size since I won't use cross fitting as usual
tau = max(data_debug[,1]) + 1  # set maximum follow-up time
#tau = 1
data_debug[,2] = data_debug[,2]*(data_debug[,1] < tau)
data_debug[,1] = ifelse(data_debug[,1] < tau, data_debug[,1], tau)
# check censoring rate 
cen.rate = mean(data_debug[,2])
print(cen.rate)
# check overlap assumption
#beta.true = c(-1, 2)  # notice the switch
beta.true = c(1, 1.5)  # notice the switch
Lambda.overlap = exp(as.matrix(data_debug[,-(1:2)])%*%beta.true) * tau
surv.overlap = exp(-Lambda.overlap)
overlap = sum(surv.overlap == 0)
print(overlap)
# split training and testing 
data_debug_train = data_debug[-10,]  # 4999 x 4 matrix 
data_debug_test = data_debug[10,]  # 4 x 1 vector 


# true survival function
X.sort = sort(unique(data_debug_train[,1]))
Lambda.true = exp(matrix(data_debug_test[-(1:2)], 1, 2)%*%beta.true) %*% t(X.sort)# baseline hazard is 1 
surv.true = exp(-Lambda.true)

# estimated survival function using Cox 
model.cox = coxph(Surv(X, Delta) ~ ., timefix = F, data = as.data.frame(data_debug_train))
Lambda.cox = exp(matrix(data_debug_test[-(1:2)], 1, 2)%*%model.cox$coefficients) %*% t(basehaz(model.cox, centered = F)$hazard)
surv.est.cox = exp(-Lambda.cox)

# plot comparison 
plot1_name <- "censor_cox_est.png"
plot1_path <- file.path(file_path, plot1_name)
png(filename = plot1_path)
plot(X.sort, surv.true, col = 1, type = "l",lty = 1,
     main = "survival curves",
     xlab = "time", ylab = "Survival probability")
lines(X.sort, surv.est.cox, col = 2, type = "l", lty = 2, lwd = 2.5)
legend("topright", 
       legend = c("True Survival", "Estimated Cox"),
       col = 1:2, lty = 1:2, lwd = c(1, 2.5), bty = "n")
dev.off()

# true survival function
X.sort.full = sort(unique(data_debug[,1]))
Lambda.true.long = exp(matrix(data_debug_test[-(1:2)], 1, 2)%*%beta.true) %*% t(X.sort.full)# baseline hazard is 1 
surv.true.long = exp(-Lambda.true.long)

# estimated survival function using Spline 
X = data_debug_train[,1]
Delta = data_debug_train[,2]
L = data_debug_train[,3]
U = data_debug_train[,4]
model.spline = hare(X, Delta, cbind(L,U))
surv.est.spline = 1 - sapply(X.sort.full, function(x) phare(q = x, cov = data_debug_test[-(1:2)], fit = model.spline))

# plot comparison 
plot2_name <- "censor_spline_est.png"
plot2_path <- file.path(file_path, plot2_name)
png(filename = plot2_path)
plot(X.sort.full, surv.true.long, col = 1, type = "l",lty = 1,
     main = "survival curves",
     xlab = "time", ylab = "Survival probability")
lines(X.sort.full, surv.est.spline, col = 3, type = "l", lty = 3, lwd = 2.5)
legend("topright", 
       legend = c("True Survival", "Estimated Spline"),
       col = c(1,3), lty = c(1,3), lwd = c(1, 2.5), bty = "n")
dev.off()

# estimated survival function using RSF 
model.rsf = rfsrc(Surv(X, Delta) ~ ., data = as.data.frame(data_debug_train), ntree = 1000, ntime = 0, mtry = 2, splitrule = 'bs.gradient')
rsf.time = model.rsf$time.interest
newdat = as.data.frame(matrix(data_debug_test, 1, 4))
colnames(newdat) = c('X', 'Delta', 'L', 'U')
surv.est.rsf = predict(model.rsf, newdata = newdat)$survival

# true survival function
Lambda.true.rsf.time = exp(matrix(data_debug_test[-(1:2)], 1, 2)%*%beta.true) %*% t(rsf.time)# baseline hazard is 1 
surv.true.rsf.time = exp(-Lambda.true.rsf.time)

# plot comparison 
plot3_name <- "censor_rsf_est.png"
plot3_path <- file.path(file_path, plot3_name)
png(filename = plot3_path)
plot(rsf.time, surv.true.rsf.time, col = 1, type = "l",lty = 1,
     main = "survival curves",
     xlab = "time", ylab = "Survival probability")
lines(rsf.time, surv.est.rsf, col = 4, type = "l", lty = 4, lwd = 2.5)
legend("topright", 
       legend = c("True Survival", "Estimated RSF"),
       col = c(1,4), lty = c(1,4),lwd = c(1, 2.5), bty = "n")
dev.off()






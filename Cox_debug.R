# this file is to check if the formula I use to estimate the nuisance with Cox 
# is correct. I use two different ways 1) calculate directly as in sim_full.R
# and 2) use predict function in Cox. These two return nearly the same plots. 
# I also use step function in R to plot, no different from regular plot().

# to ensure the prediction is accurate, given the data generating process is Cox, 
# I also draw the true survival function at the end and compare. 
# since the true baseline hazard model is lambda0(t) = 1, Lambda(t) = int_0^t 1 dv = t
# and the coefficient are beta.true = c(-1,2). The result shows that the Cox fomula and 
# predict function in Cox are well performed. 

# in summary, since the Cox formula implemented in Jiyu perform correctly, I will assume
# his spline and RSF are also correct. 

set.seed(2007)
n = 1000
n.original=n
n = ceiling(n*1.01)
eps_1 = runif(n,-1,1)
eps_2 = runif(n)
eps_3 = runif(n)
L = 0.5*eps_1 + rnorm(n,0, sd = 1)  # mean zero 
U = eps_1^2 + rnorm(n,0, sd = 0.3)  # mean 1/3
Z = cbind(L,U)
t = -log(eps_2)/exp(as.vector(Z%*%c(-1, 2)))
c = -log(eps_3)/exp(as.vector(Z%*%c(0.5, 2.5)))
X = pmin(t,c)
Delta = as.numeric(t < c)
data = cbind(t, c, X, Delta, L, U)[1:n.original,]
data = data[,-c(1,2)]
tau = 3 # set maximum follow-up time 
data[,2] = data[,2]*(data[,1] < tau)
data[,1] = ifelse(data[,1] < tau, data[,1], tau)

data = as.data.frame(data)
X.full = data[,1]
X.sort.full = sort(unique(X.full))
folds = cut(1:nrow(data), breaks = 2, labels = F)
S_t.temp = S_c.temp = matrix(0, nrow = nrow(data), ncol = length(unique(X.full)))

fold=1
dfs = dfs.c = data[folds==fold,]
dfn = dfn.c = data[folds!=fold,]
X = dfn[,1]
X.sort = sort(unique(X))
Delta = dfn[,2]
Delta_c = (1 - Delta)*(X < tau)
dfn.c$Delta = Delta_c
colnames(dfn.c)[2] = 'Delta_c'
colnames(dfs.c)[2] = 'Delta_c'
L = dfn[,3]
U = dfn[,4]

model = coxph(Surv(X, Delta) ~ ., timefix = F, data = dfn)
print(model$coefficients)

model_c = coxph(Surv(X, Delta_c) ~ ., timefix = F, data = dfn.c)
print(model_c$coefficients)

# predict using formula 
basehaz = basehaz(model, centered = F)
Lambda = exp(as.matrix(dfs[,-(1:2)])%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
surv = exp(-Lambda)

# predict using cox predict function 
newdat = data.frame(X = basehaz$time,
                    Delta = 1,
                    L = rep(dfs[1, 3], length(basehaz$time)),
                    U = rep(dfs[1, 4], length(basehaz$time)))
pred <- predict(model, newdata = newdat, type = "survival", se.fit = TRUE)
pred_surv = pred$fit

# create two step functions 
stepf = stepfun(basehaz$time, c(1, surv[1,]))
pred_surv_stepf = stepfun(newdat$X, c(1,pred_surv))

# calculate true survival probability 
beta.true = c(-1, 2)
Lambda.true = exp(as.matrix(dfs[,-(1:2)])%*%beta.true) %*% t(basehaz$time)# baseline hazard is 1 
surv.true = exp(-Lambda.true)

# plot of predicted surv through predict function in Cox 
plot(pred_surv_stepf, do.points = FALSE, col = 1, lty = 1,
     main = "Predicted survival curves",
     xlab = "Days", ylab = "Survival probability")
# plot of predicted surv through direction formula  
plot(stepf, do.points = FALSE, col = 2, lty = 1, add = TRUE)
# plot of predicted surv through direction formula  
lines(basehaz$time, surv[1,], col = 3, type = "l")
lines(basehaz$time, surv.true[1,], col = 4, type = "l")






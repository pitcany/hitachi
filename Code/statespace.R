library(readr)
cmapps <- read_table2("Documents/Research/bayesian_control/CMAPPS/train_FD001.txt",col_names = FALSE)
colnames(cmapps)[1] <- 'machine'
colnames(cmapps)[2] <- 'timestep'

# Load Libraries ---- 
library(dlm)
#library(tidyverse) # for plotting

library(astsa)
# Setup
y = cbind(globtemp, globtempl); num = nrow(y); input = rep(1,num)
A = array(rep(1,2), dim=c(2,1,num))
mu0 = -.35; Sigma0 = 1; Phi = 1

# Function to evaluate the likelihood
Linn = function(para){
  cQ = para[1] # sigma_w
  cR1 = para[2] # 11 element of chol(R)
  cR2 = para[3] # 22 element of chol(R)
  cR12 = para[4] # 12 element of chol(R)
  cR = matrix(c(cR1,0,cR12,cR2),2) # build matrix
  drift = para[5]
  kf = Kfilter1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
  return (kf$like) }
  #Estimation
init.par = c(.1,.1,.1,0,.05) # initial values of parameters
(est = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE,
             control=list(trace=1,REPORT=1))) # output not shown
SE = sqrt(diag(solve(est$hessian)))
# Display estimates
u = cbind(estimate=est$par, SE)
rownames(u)=c('sigw','cR11','cR22','cR12','drift'); u
# Smooth (first set parameters to their final estimates)
cQ = est$par[1]
cR1 = est$par[2]
cR2 = est$par[3]
cR12 = est$par[4]
cR = matrix(c(cR1,0,cR12,cR2), 2)
(R = t(cR)%*%cR) # to view the estimated R matrix
drift = est$par[5]
ks = Ksmooth1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
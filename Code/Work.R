library(readr)
cmapps <- read_table2("Documents/Research/bayesian_control/CMAPPS/train_FD001.txt",col_names = FALSE)
colnames(cmapps)[1] <- 'machine'
colnames(cmapps)[2] <- 'timestep'

# Load Libraries ---- 
# library(dlm)
#library(tidyverse) # for plotting

library(astsa)
# Setup

y <- cmapps[cmapps$machine==1,6:26]
y[c("X6","X10","X11","X13","X15","X18","X21","X23","X24")]<-NULL
#num <- nrow(y)
num <- nrow(y)
A=array(1, dim=c(12,4))
Phi = array(rep(.5,16),dim=c(4,4))
(cR <- matrix(rep(1,144), 12, 12))
# lower.tri(cR)
cR[lower.tri(cR)] <- 0
cQ = matrix(c(1,0,0,0,1,1,0,0,1,1,1,0,1,1,1,1),nrow=4,ncol=4)
Sigma0= array(rep(.5,16),dim=c(4,4))
mu0=c(1,0,0,0)

# set.seed(999); num = 100
# x = arima.sim(n=num+1, list(ar=.8), sd=1)
# y = ts(x[-1] + rnorm(num,0,1))
# u = ts.intersect(y, lag(y,-1), lag(y,-2))
# varu = var(u); coru = cor(u)
# phi = coru[1,3]/coru[1,2]
# q = (1-phi^2)*varu[1,2]/phi
# r = varu[1,1] - q/(1-phi^2)

# naive initialization -- vector quantization

#(em = EM0(num, y, A=1, mu0=0, Sigma0=2.8, Phi=phi, cQ=sqrt(q), cR=sqrt(r),
#          max.iter=75, tol=.00001))
#astsa::Kfilter0(num, y, A, mu0, Sigma0, Phi, cQ, cR) 

# EM0_altered(num, y, A, mu0, Sigma0, Phi, cQ, cR,
#             max.iter=75, tol=.00001)
#(em = EM0_altered(num, y, A, mu0, Sigma0, Phi, cQ, cR,
#         max.iter=75, tol=.00001))

# EM procedure - output not show
EM0_altered <- function( num,y,A,mu0,Sigma0,Phi,cQ,cR,max.iter=50,tol=.01 ){
    Phi=as.matrix(Phi)
    pdim=nrow(Phi)
    y=as.matrix(y)
    qdim=ncol(y)
    cvg=1+tol
    like=matrix(0,max.iter,1)
    cat("iteration","   -loglikelihood", "\n")
    #----------------- start EM -------------------------
    for(iter in 1:max.iter){ 
      ks=astsa::Ksmooth0(num,y,A,mu0,Sigma0,Phi,cQ,cR)
      like[iter]=ks$like
      cat("   ",iter, "        ", ks$like, "\n")
      if(iter>1) cvg=(like[iter-1]-like[iter])/abs(like[iter-1])
      if(cvg<0) stop("Likelihood Not Increasing")
      if(abs(cvg)<tol) break
      # Lag-One Covariance Smoothers 
      Pcs=array(NA, dim=c(pdim,pdim,num))     # Pcs=P_{t,t-1}^n
      eye=diag(1,pdim)
      Pcs[,,num]=(eye-ks$Kn%*%A)%*%Phi%*%ks$Pf[,,num-1] #
      for(k in num:3){
        Pcs[,,k-1]=ks$Pf[,,k-1]%*%t(ks$J[,,k-2])+
          ks$J[,,k-1]%*%(Pcs[,,k]-Phi%*%ks$Pf[,,k-1])%*%t(ks$J[,,k-2])}
      Pcs[,,1]=ks$Pf[,,1]%*%t(ks$J0)+
        ks$J[,,1]%*%(Pcs[,,2]-Phi%*%ks$Pf[,,1])%*%t(ks$J0)       
      # Estimation
      S11 = ks$xs[,,1]%*%t(ks$xs[,,1]) + ks$Ps[,,1]
      S10 = ks$xs[,,1]%*%t(ks$x0n) + Pcs[,,1]
      S00 = ks$x0n%*%t(ks$x0n) + ks$P0n
      u = y[1,]-A%*%ks$xs[,,1]
      R = u%*%t(u) + A%*%ks$Ps[,,1]%*%t(A)
      # for observation matrix update
      E = y[1,]%*%t(ks$xs[,,1])
      for(i in 2:num){
        S11 = S11 + ks$xs[,,i]%*%t(ks$xs[,,i]) + ks$Ps[,,i]
        S10 = S10 + ks$xs[,,i]%*%t(ks$xs[,,i-1]) + Pcs[,,i]
        S00 = S00 + ks$xs[,,i-1]%*%t(ks$xs[,,i-1]) + ks$Ps[,,i-1]
        u = y[i,]-A%*%ks$xs[,,i]
        E = E + y[i,]%*%t(ks$xs[,,i])
        R = R + u%*%t(u) + A%*%ks$Ps[,,i]%*%t(A)
        svd(R)
      }
      # A is observation matrix updated
      A = E %*% solve(S11+ks$x0n%*%t(ks$x0n) + ks$P0n)
      Phi=S10%*%solve(S00)
      Q=(S11-Phi%*%t(S10))/num
      Q=(t(Q)+Q)/2        # make sure symmetric
      cQ=chol(Q)
      R=R/num
      cR=chol(R)
      mu0=ks$x0n
      # mu0=mu0              # uncomment this line to keep mu0 fixed
      Sigma0=ks$P0n
      # Sigma0=Sigma0        # uncomment this line to keep Sigma0 fixed
    }
    list(A=A,Phi=Phi,Q=Q,R=R,mu0=mu0,Sigma0=Sigma0,like=like[1:iter],niter=iter,cvg=cvg)
  }

em <- EM0_altered(num, y, A, mu0, Sigma0, Phi, cQ, cR,
            max.iter=75, tol=.00001)
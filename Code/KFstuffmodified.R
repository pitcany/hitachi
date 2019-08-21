
Kfilter0 <-
  function(num,y,A,mu0,Sigma0,Phi,cQ,cR){
    #
    # NOTE: must give cholesky decomp: cQ=chol(Q), cR=chol(R)
    Q=t(cQ)%*%cQ
    R=t(cR)%*%cR
    # y is num by q  (time=row series=col)
    # A is a q by p matrix
    # R is q by q
    # mu0 is p by 1
    # Sigma0, Phi, Q are p by p
    Phi=as.matrix(Phi)
    pdim=nrow(Phi)    
    y=as.matrix(y)
    qdim=ncol(y)
    xp=array(NA, dim=c(pdim,1,num))         # xp=x_t^{t-1}          
    Pp=array(NA, dim=c(pdim,pdim,num))      # Pp=P_t^{t-1}
    xf=array(NA, dim=c(pdim,1,num))         # xf=x_t^t
    Pf=array(NA, dim=c(pdim,pdim,num))      # Pf=x_t^t
    innov=array(NA, dim=c(qdim,1,num))      # innovations
    sig=array(NA, dim=c(qdim,qdim,num))     # innov var-cov matrix
    # initialize (because R can't count from zero)
    x00=as.matrix(mu0, nrow=pdim, ncol=1)
    P00=as.matrix(Sigma0, nrow=pdim, ncol=pdim)
    xp[,,1]=Phi%*%x00
    Pp[,,1]=Phi%*%P00%*%t(Phi)+Q
    sigtemp=A%*%Pp[,,1]%*%t(A)+R
    sig[,,1]=(t(sigtemp)+sigtemp)/2     # innov var - make sure it's symmetric
    siginv=solve(sig[,,1])          
    K=Pp[,,1]%*%t(A)%*%siginv
    innov[,,1]=y[1,]-A%*%xp[,,1]
    xf[,,1]=xp[,,1]+K%*%innov[,,1]
    Pf[,,1]=Pp[,,1]-K%*%A%*%Pp[,,1]
    sigmat=as.matrix(sig[,,1], nrow=qdim, ncol=qdim)
    like = log(det(sigmat)) + t(innov[,,1])%*%siginv%*%innov[,,1]   # -log(likelihood)
    ########## start filter iterations ###################
    for (i in 2:num){
      if (num < 2) break
      xp[,,i]=Phi%*%xf[,,i-1]
      Pp[,,i]=Phi%*%Pf[,,i-1]%*%t(Phi)+Q
      sigtemp=A%*%Pp[,,i]%*%t(A)+R
      sig[,,i]=(t(sigtemp)+sigtemp)/2     # innov var - make sure it's symmetric
      siginv=solve(sig[,,i])              
      K=Pp[,,i]%*%t(A)%*%siginv
      innov[,,i]=y[i,]-A%*%xp[,,i]
      xf[,,i]=xp[,,i]+K%*%innov[,,i]
      Pf[,,i]=Pp[,,i]-K%*%A%*%Pp[,,i]
      sigmat=as.matrix(sig[,,i], nrow=qdim, ncol=qdim)
      like= like + log(det(sigmat)) + t(innov[,,i])%*%siginv%*%innov[,,i]
    }
    like=0.5*like
    list(xp=xp,Pp=Pp,xf=xf,Pf=Pf,like=like,innov=innov,sig=sig,Kn=K)
  }

Ksmooth0 <-
  function(num,y,A,mu0,Sigma0,Phi,cQ,cR){
    #
    # Note: Q and R are given as Cholesky decomps
    #       cQ=chol(Q), cR=chol(R)
    #
    kf=astsa::Kfilter0(num,y,A,mu0,Sigma0,Phi,cQ,cR)
    pdim=nrow(as.matrix(Phi))  
    xs=array(NA, dim=c(pdim,1,num))      # xs=x_t^n
    Ps=array(NA, dim=c(pdim,pdim,num))   # Ps=P_t^n
    J=array(NA, dim=c(pdim,pdim,num))    # J=J_t
    xs[,,num]=kf$xf[,,num] 
    Ps[,,num]=kf$Pf[,,num]
    for(k in num:2)  {
      J[,,k-1]=(kf$Pf[,,k-1]%*%t(Phi))%*%solve(kf$Pp[,,k])
      xs[,,k-1]=kf$xf[,,k-1]+J[,,k-1]%*%(xs[,,k]-kf$xp[,,k])
      Ps[,,k-1]=kf$Pf[,,k-1]+J[,,k-1]%*%(Ps[,,k]-kf$Pp[,,k])%*%t(J[,,k-1])
    }
    # and now for the initial values because R can't count backward to zero
    x00=mu0
    P00=Sigma0
    J0=as.matrix((P00%*%t(Phi))%*%solve(kf$Pp[,,1]), nrow=pdim, ncol=pdim)
    x0n=as.matrix(x00+J0%*%(xs[,,1]-kf$xp[,,1]), nrow=pdim, ncol=1)
    P0n= P00 + J0%*%(Ps[,,1]-kf$Pp[,,1])%*%t(J0)
    list(xs=xs,Ps=Ps,x0n=x0n,P0n=P0n,J0=J0,J=J,xp=kf$xp,Pp=kf$Pp,xf=kf$xf,Pf=kf$Pf,like=kf$like,Kn=kf$K)
  }

EM0 <-
  function(num,y,A,mu0,Sigma0,Phi,cQ,cR,max.iter=50,tol=.01){
    #
    # Note: Q and R are given as Cholesky decomps
    #       cQ=chol(Q), cR=chol(R)
    #
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
      Pcs[,,num]=(eye-ks$Kn%*%A)%*%Phi%*%ks$Pf[,,num-1]
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
      for(i in 2:num){
        S11 = S11 + ks$xs[,,i]%*%t(ks$xs[,,i]) + ks$Ps[,,i]
        S10 = S10 + ks$xs[,,i]%*%t(ks$xs[,,i-1]) + Pcs[,,i]
        S00 = S00 + ks$xs[,,i-1]%*%t(ks$xs[,,i-1]) + ks$Ps[,,i-1]
        u = y[i,]-A%*%ks$xs[,,i]
        R = R + u%*%t(u) + A%*%ks$Ps[,,i]%*%t(A)
      }
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
    list(Phi=Phi,Q=Q,R=R,mu0=mu0,Sigma0=Sigma0,like=like[1:iter],niter=iter,cvg=cvg)
  }


library(matrixStats)
library(parallel)
library(Rcpp)
library(abind)
library(scales)
Rcpp::sourceCpp('MIS_V3.cpp')
source('MIS_V3_target.R')

logUnPi<-function(x){
  
 lDen<-logTargetUnDens(x,lUnDens_list,lWeights)
  return(lDen)
}

rPi<-function(N,d){
  
  mu_vec_Pi<-c(-5,0,10)
  sd_vec_Pi<-c(2,5,0.5)
  indDist<-sample(1:3,N,rep=T)
  
  mu_sample<-mu_vec_Pi[indDist]
  sd_sample<-sd_vec_Pi[indDist]
  
  out<-sapply(1:N, function(n) rnorm(d,mean=mu_sample,sd=sd_sample) )
  return( out )
}

logGamma<-function(x,mu_ln,sigma_ln,mu_vec_sub,sigma_vec_sub,N_vec,beta_t){
  
  #logNum<-beta_t*logUnPi(x)+logQ_L(x,mu_ln)
  logNum<-beta_t*logUnPi(x)+logQ_L_cpp(x,mu_ln,sigma_ln)
  
  #tmp<-sapply(1:nrow(mu_vec_sub),function(i) beta_t*logQ_L(x,mu_vec_sub[i,]) + log(N_vec[i]) )
  #tmp<-logGamma_den_cpp(x,mu_vec_sub,N_vec,rep(sigma_prop,d),beta_t)
  logQ_vec<-logQ_L_vec_cpp(x,mu_vec_sub,sigma_vec_sub)
  tmp<-beta_t*logQ_vec+log(N_vec)
  
  logDen<-logSumExp(tmp)
  
  return(logNum-logDen) 
}

logQ_L<-function(x,mu_l,sigma_l){
  
  d<-length(x)
  #out<-sum(dnorm(x,mean=mu_l,sd=sigma_prop,log=T))
  out<-logQ_L_cpp(x,mu_l,sigma_l)
  return(out)
}

randQ_L<-function(mu_l,sigma_l){
  
  mu<-as.double(mu_l)
  
  out<-rnorm(d,mean=mu,sd=sigma_l)
  return(out)
}

Kernel_n<-function(x,mu_ln,sigma_ln,mu_vec,sigma_vec,N_vec,t){
  
  d<-length(x)
  x_star<-x+rnorm(d,0,sigma_MCMC)
  lratio<-logGamma(x_star,mu_ln,sigma_ln,mu_vec,sigma_vec,N_vec,t)-logGamma(x,mu_ln,sigma_ln,mu_vec,sigma_vec,N_vec,t)
  
  if(log(runif(1))<=lratio){
    
    accepted<-1
    return( c(accepted,x_star) )
  }
  else{
    
    accepted<-0
    return( c(accepted,x) )
  }
}

Kernel<-function(x,mu_vec,sigma_vec,N_vec,t){
  
  sapply(1:N,)
  
  d<-length(x)
  x_star<-x+rnorm(d,0,sigma_MCMC)
  lratio<-logGamma(x_star,mu_ln,sigma_ln,mu_vec,sigma_vec,N_vec,t)-logGamma(x,mu_ln,sigma_ln,mu_vec,sigma_vec,N_vec,t)
  
  if(log(runif(1))<=lratio){
    
    accepted<-1
    return( c(accepted,x_star) )
  }
  else{
    
    accepted<-0
    return( c(accepted,x) )
  }
}

system_resamp<-function(N,logW){
  
  U<-runif(1,0,1/N)
  U_i<-((1:N)-1)/N+U
  
  w<-exp(logW-max(logW))
  w_cum<-cumsum(w)/sum(w)
  
  indices<-sapply(U_i,function(u){
    
    out<-min(which(w_cum>=u))
    return(out)
  })
  
  return(indices)
}

logESS<-function(lW){
  
  lWnorm<-lW-logSumExp(lW)
  
  return(-logSumExp(2*lWnorm))
}

log_relCESS<-function(lW,lWinc){
  
  lWnorm<-lW-logSumExp(lW)
  lNum<-2*logSumExp(lWnorm+lWinc)
  lDen<-logSumExp(lWnorm+2*lWinc)
  
  return(lNum-lDen)
}


MIS<-function(d,mu_vec,sigma_vec,N,plots=F){
  
  beta_adap<-c(0)
  
  L<-matrix(0,nrow=N,ncol=1)
  L[,1]<-sample(1:nrow(mu_vec),N,rep=T)
  
  tab<-table(L[,1])
  ind_vec0<-as.double(names(tab))
  if(d==1){
    mu_vec_sub<-matrix(mu_vec[ind_vec0,],ncol=1)
    sigma_vec_sub<-matrix(sigma_vec[ind_vec0,],ncol=1)
  }
  else{
    mu_vec_sub<-mu_vec[ind_vec0,]
    sigma_vec_sub<-sigma_vec[ind_vec0,]
  }
  N_vec0<-as.double(tab)
  
  X<-array(0,dim=c(d,N,1))
  X[,,1]<-sapply(1:N,function(n) randQ_L(mu_vec[L[n,1],],sigma_vec[L[n,1]]))
  if(plots){
    # plot(density(rPi(5e2)),col="blue",lwd=2,xlim=c(min(mu_vec)*(1-sign(min(mu_vec))*.2),max(mu_vec)*(1+sign(max(mu_vec))*.2)))
    # lines(density(X[,1]),col="red",lwd=2)
  }
  
  logW_BH<-sapply(1:N,function(n) logGamma(X[,n,1],mu_vec[L[n,1],],sigma_vec[L[n,1],],mu_vec_sub,sigma_vec_sub,N_vec0,1) - 
                      logGamma(X[,n,1],mu_vec[L[n,1],],sigma_vec[L[n,1],],mu_vec_sub,sigma_vec_sub,N_vec0,0) )
  
  logW_RB<-sapply(1:N,function(n) logUnPi(X[,n,1]) - 
                       logSumExp( sapply(1:nrow(mu_vec),function(l) logQ_L_cpp(X[,n,1],mu_vec[l,],sigma_vec[l,])) -log(nrow(mu_vec)) ))
  
  logZ_BH<-logSumExp(logW_BH)-log(N)
  logZ_RB<-logSumExp(logW_RB)-log(N)
  
  return(list(logW_BH=logW_BH,logW_RB=logW_RB,X=X[,,1],L=L,logZ_BH=logZ_BH,logZ_RB=logZ_RB))
}

sAIS_MIS<-function(d,mu_vec,sigma_vec,N,KernSteps,beta_thresh,plots=F,beta_vec){
  
  beta_adap<-beta_vec
  
  L<-matrix(0,nrow=N,ncol=1)
  L[,1]<-sample(1:nrow(mu_vec),N,rep=T)
  #L[,1]<-sample(1:K,K,rep=F)
  tab<-table(L[,1])
  ind_vec0<-as.double(names(tab))
  if(d==1){
    mu_vec_sub<-matrix(mu_vec[ind_vec0,],ncol=1)
    sigma_vec_sub<-matrix(sigma_vec[ind_vec0,],ncol=1)
  }
  else{
    mu_vec_sub<-mu_vec[ind_vec0,]
    sigma_vec_sub<-sigma_vec[ind_vec0,]
  }
  N_vec0<-as.double(tab)
  
  X<-array(0,dim=c(d,N,2))
  X[,,1]<-sapply(1:N,function(n) randQ_L(mu_vec[L[n,1],],sigma_vec[L[n,1]]))
  if(plots){
    # plot(density(rPi(5e2)),col="blue",lwd=2,xlim=c(min(mu_vec)*(1-sign(min(mu_vec))*.2),max(mu_vec)*(1+sign(max(mu_vec))*.2)))
    # lines(density(X[,1]),col="red",lwd=2)
  }
  
  logW<-matrix(0,nrow=N,ncol=1)
  logW[,1]<-0
  logWinc<-matrix(0,nrow=N,ncol=0)
  A<-matrix(0,nrow=N,ncol=0)
  
  logZ<-0
  #for(t in 1:(nDists)){
  t<-1
  while(beta_adap[t]<1){
    
    logWinc<-cbind(logWinc,rep(0,N))
    logW<-cbind(logW,rep(0,N))
    A<-cbind(A,rep(0,N))
    L<-cbind(L,rep(0,N))
    X<-abind(X, array(0, replace(dim(X), 3, t+2)), along = 3)
    
    logG_t<-sapply(1:N,function(n) logGamma(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t]) )
    logQ_vec_den<-sapply(1:N,function(n) logQ_L_vec_cpp(X[,n,t],mu_vec_sub,sigma_vec_sub) )
    logQ_vec_num<-sapply(1:N,function(n) logQ_L_cpp(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],]) )
    logUnPi_num<-sapply(1:N,function(n) logUnPi(X[,n,t]) )
    
    logWinc_t<-sapply(1:N,function(n) logGamma(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t+1]) ) -
      logSumExp(logG_t)
    
    logWinc[,t]<-logWinc_t
    
    logW[,t+1]<-logW[,t]+logWinc_t
    
    l_rCESS_t<-log_relCESS(logW[,t],logWinc_t)
    if(plots)
      print(paste0("beta: ",beta_adap[t+1],", relCESS: ",exp(l_rCESS_t)))
    
    lESS_t<-logESS(logW[,t+1])
    if(plots)
      print(paste0("relESS: ",exp(lESS_t-log(N))))
    
    A_tp1<-1:N
    
    X_tp1<-X[,A_tp1,t]
    if(d==1)
      X_tp1<-matrix(X_tp1,nrow=1)
    L_tp1<-L[A_tp1,t]
    
    logZ<-logZ+logSumExp(logW[,t+1])
    logW[,t+1]<-0
    
    if(KernSteps>0){
      
      for(i in 1:KernSteps){

        log_pn<-mapply(function(n) logGamma(X_tp1[,n],mu_vec[L_tp1[n],],sigma_vec[L_tp1[n],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t+1]) -
                         logQ_L_cpp(X_tp1[,n],mu_vec[L_tp1[n],],sigma_vec[L_tp1[n],]) ,1:N)
        n_star<-sample(N,1,prob = exp(log_pn-max(log_pn)))
        X_tp1[,-n_star]<-sapply(1:N,function(n) randQ_L(mu_vec[L_tp1[n],],sigma_vec[L_tp1[n],]))[,-n_star]

      }
    }
    
    # if(plots)
    #   lines(density(X_tp1),lty=2)
    
    A[,t]<-A_tp1
    X[,,t+1]<-X_tp1
    L[,t+1]<-L_tp1
    
    t<-t+1
  }
  
  return(list(logW=logW[,t],X=X,L=L,A=A,logZ=logZ,beta_vec=beta_adap))
}

AIS_MIS<-function(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh,beta_thresh,plots=F,beta_vec=NULL){
  
  if(is.null(beta_vec)){
    beta_adap<-c(0)
    adapt<-1
  }
  else{
    beta_adap<-beta_vec
    adapt<-0
  }
  
  L<-matrix(0,nrow=N,ncol=1)
  L[,1]<-sample(1:nrow(mu_vec),N,rep=T)
  #L[,1]<-sample(1:K,K,rep=F)
  tab<-table(L[,1])
  ind_vec0<-as.double(names(tab))
  if(d==1){
    mu_vec_sub<-matrix(mu_vec[ind_vec0,],ncol=1)
    sigma_vec_sub<-matrix(sigma_vec[ind_vec0,],ncol=1)
  }
  else{
    mu_vec_sub<-mu_vec[ind_vec0,]
    sigma_vec_sub<-sigma_vec[ind_vec0,]
  }
  N_vec0<-as.double(tab)
  
  X<-array(0,dim=c(d,N,2))
  X[,,1]<-sapply(1:N,function(n) randQ_L(mu_vec[L[n,1],],sigma_vec[L[n,1]]))
  if(plots){
    # plot(density(rPi(5e2)),col="blue",lwd=2,xlim=c(min(mu_vec)*(1-sign(min(mu_vec))*.2),max(mu_vec)*(1+sign(max(mu_vec))*.2)))
    # lines(density(X[,1]),col="red",lwd=2)
  }
  
  logW<-matrix(0,nrow=N,ncol=1)
  logW[,1]<--log(N)
  logWinc<-matrix(0,nrow=N,ncol=0)
  A<-matrix(0,nrow=N,ncol=0)
  
  logZ<-0
  #for(t in 1:(nDists)){
  t<-1
  while(beta_adap[t]<1){
    
    logWinc<-cbind(logWinc,rep(0,N))
    logW<-cbind(logW,rep(0,N))
    A<-cbind(A,rep(0,N))
    L<-cbind(L,rep(0,N))
    X<-abind(X, array(0, replace(dim(X), 3, t+2)), along = 3)
    
    logG_t<-sapply(1:N,function(n) logGamma(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t]) )
    logQ_vec_den<-sapply(1:N,function(n) logQ_L_vec_cpp(X[,n,t],mu_vec_sub,sigma_vec_sub) )
    logQ_vec_num<-sapply(1:N,function(n) logQ_L_cpp(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],]) )
    logUnPi_num<-sapply(1:N,function(n) logUnPi(X[,n,t]) )
    
    if(adapt==1){
      
      log_relCESS_fn<-function(beta){
        
        logNum<-beta*logUnPi_num+logQ_vec_num
        tmp<-beta*logQ_vec_den+log(N_vec0)
        logDen<-sapply(1:N,function(n) logSumExp(tmp[,n]) )
        
        lWinc<-logNum-logDen - logG_t
        
        out<-log_relCESS(logW[,t],lWinc)
        
        return(out-log(beta_thresh))
      }
      
      #browser()
      if(log_relCESS_fn(1)>0)
        beta_adap[t+1]<-1
      else
        beta_adap[t+1]<-uniroot(log_relCESS_fn,c(beta_adap[t],1))$root
    }
    
    logWinc_t<-sapply(1:N,function(n) logGamma(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t+1]) - 
                        logGamma(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t]) )
    logWinc[,t]<-logWinc_t
    
    logW[,t+1]<-logW[,t]+logWinc_t
    
    l_rCESS_t<-log_relCESS(logW[,t],logWinc_t)
    if(plots)
      print(paste0("beta: ",beta_adap[t+1],", relCESS: ",exp(l_rCESS_t)))
    
    lESS_t<-logESS(logW[,t+1])
    if(plots)
      print(paste0("relESS: ",exp(lESS_t-log(N))))
    
    #if(resamp_vec[t]==1 & t%%resamp_every==0){
    if(lESS_t<log(alpha_thresh)+log(N)){
      
      A_tp1<-system_resamp(N,logW[,t+1])
      
      if(plots)
        print("Resampled!")
      
      X_tp1<-X[,A_tp1,t]
      if(d==1)
        X_tp1<-matrix(X_tp1,nrow=1)
      L_tp1<-L[A_tp1,t]
      #L_tp1<-L[1:N,t]
      
      logZ<-logZ+logSumExp(logW[,t+1])
      
      logW[,t+1]<--log(N)
      
      #logW[,t+1]<-sapply(1:N,function(n) logGamma(X_tp1[n],L_tp1[n],ind_vec,N_vec,t+1) )-
      #  logSumExp( sapply(1:N, function(n) logGamma(X_tp1[n],L_tp1[n],ind_vec,N_vec,t+1) ))
    }
    else{
      
      A_tp1<-1:N
      
      X_tp1<-X[,A_tp1,t]
      if(d==1)
        X_tp1<-matrix(X_tp1,nrow=1)
      L_tp1<-L[A_tp1,t]
    }
    
    if(KernSteps>0){
      
      for(i in 1:KernSteps){
        
        tmp<-mapply(function(n) Kernel_n(X_tp1[,n],mu_vec[L_tp1[n],],sigma_vec[L_tp1[n],],mu_vec_sub,sigma_vec_sub,N_vec0,beta_adap[t+1]),1:N)
        X_tp1<-tmp[-1,]
        
        if(plots)
          print(mean(tmp[1,]))
      }
    }
    
    # if(plots)
    #   lines(density(X_tp1),lty=2)
    
    A[,t]<-A_tp1
    X[,,t+1]<-X_tp1
    L[,t+1]<-L_tp1
    
    t<-t+1
  }
  
  # if(plots)
  #   lines(density(X[,nDists+1]),lwd=2)
  
  logZ<-logZ+logSumExp(logW[,t])
  
  return(list(logW=logW[,t],X=X,L=L,A=A,logZ=logZ,beta_vec=beta_adap))
}


detectCores()
options(mc.cores=detectCores()-2)
getOption("mc.cores", 2L)

d<-5; N<-2000; KernSteps<-1; alpha_thresh<-0.5; beta_thresh<-0.99; limInf<--30; limSup<-30
K<-max(floor(N^(1/d)),2)
K<-max(floor((N/2)^(1/d)),2)
mu<-seq(limInf,limSup,length.out = 2*K+1)[2*(1:K)]
sigma_MCMC<-max(diff(mu)/2)/sqrt(d)
mu_vec<-as.matrix(expand.grid(replicate(d,mu,simplify = F)))
sigma_vec<-matrix(rep(max(diff(mu)/2),d),nrow=nrow(mu_vec),ncol=d)

nProps<-50
if(nProps==1){
  x1<-0
  x2<-.03*(100-x1^2)
}else {
  x1<-seq(limInf,limSup,length.out=nProps)
  x2<-.03*(100-x1^2)
}
mu_vec<-matrix(0,nrow=50,ncol=d)
sigma_vec<-matrix(0,nrow=50,ncol=d)
for(i in 1:nProps){
  
  mu_vec[i,]<-c(x1[i],x2[i],rep(0,d-2))
  sigma_vec[i,]<-c(2,2,rep(2,d-2))
}

source('MIS_V3_target.R')
lUnDens_list<-list(ld0,ld1,ld2,ld3)
lWeights<-c(-Inf,-0,-Inf,-Inf)
lZ_true_vec<-c(ld0(0,"lZ",d),ld1(0,"lZ",d),ld2(0,"lZ",d),ld3(0,"lZ",d))
lZ_true<-logSumExp(lZ_true_vec+lWeights)-logSumExp(lWeights)

MIS_results<-MIS(d,mu_vec,sigma_vec,N,plots=T)

indRes_BH<-sample(N,rep=T,prob=exp(MIS_results$logW_BH-max(MIS_results$logW_BH)))
indRes_RB<-sample(N,rep=T,prob=exp(MIS_results$logW_RB-max(MIS_results$logW_RB)))
plot(mu_vec,xlab = "x1", ylab = "x2")
points(mu_vec[MIS_results$L[indRes_BH,ncol(MIS_results$L)],],pch=16,col=alpha("black",0.25))

i1<-1; i2<-2
x1<-seq(limInf,limSup,.1); x2<-x1
targetUnDens<-function(x,y){
  
  theta<-rep(0,d)
  theta[c(i1,i2)]<-c(x,y)
  
  return(logTargetUnDens(theta,lUnDens_list,lWeights))
}
z<-outer(x1,x1,FUN=Vectorize(targetUnDens))
contour(x1,x2, z, drawlabels=FALSE, nlevels=100, add=T)

plot(t(MIS_results$X[,indRes_BH]), xlab = "x1", ylab = "x2")
plot(t(MIS_results$X[,indRes_RB]), xlab = "x1", ylab = "x2")
print(c(lZ_true,MIS_results$logZ_RB,MIS_results$logZ_BH))


AIS_results<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.9,T)
print(c(lZ_true,AIS_results$logZ))
beta_vec<-AIS_results$beta_vec

plot(matrix(1:ncol(AIS_results$L),nrow=N,ncol=ncol(AIS_results$L),byrow = T),AIS_results$L)

indRes_AIS<-sample(N,rep=T,prob=exp(AIS_results$logW-max(AIS_results$logW)))
plot(mu_vec,xlab = "x1", ylab = "x2")
points(mu_vec[AIS_results$L[indRes_AIS,ncol(AIS_results$L)],],pch=16,col=alpha("black",0.25))
contour(x1,x2, z, drawlabels=FALSE, nlevels=100, add=T)
plot(t(AIS_results$X[,indRes_AIS,ncol(AIS_results$L)]))

sAIS_results<-sAIS_MIS(d,mu_vec,sigma_vec,N,KernSteps=0,beta_thresh,T,AIS_results$beta_vec)
print(c(lZ_true,sAIS_results$logZ))

rep<-100
logZ_vec<-mcmapply(function(i){
  
  nDists<-1
  tmp_MIS<-MIS(d,mu_vec,sigma_vec,N)
  
  return(c(tmp_MIS$logZ_RB,tmp_MIS$logZ_BH))
},1:rep)
boxplot(t(logZ_vec),names=c("RB","BH"),ylab="logZ"); abline(h=lZ_true)

plot(density(logZ_vec[2,]),lty=2)
lines(density(logZ_vec[1,]))
sapply(1:nrow(logZ_vec),function(n) var(exp(logZ_vec[n,])))

rep<-100
logZ_vec<-mcmapply(function(i){
  
  tmp_MIS<-MIS(d,mu_vec,sigma_vec,N)
  KernSteps<-1;
  tmp_BH_sAIS<-sAIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,beta_thresh=0.9,T,beta_vec)
  tmp_BH_AIS<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0,beta_thresh=0.9,beta_vec)
  tmp_BH_Boot<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=1,beta_thresh=0.9,beta_vec)
  tmp_BH_SMC<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.9,beta_vec)
  
  return(c(tmp_MIS$logZ_RB,tmp_MIS$logZ_BH,tmp_BH_sAIS$logZ,tmp_BH_AIS$logZ,tmp_BH_Boot$logZ,tmp_BH_SMC$logZ))
},1:rep)

length(beta_vec)
boxplot(t(logZ_vec),names=c("RB","BH","sAIS","mAIS","Boot","Adap")); abline(h=lZ_true)
boxplot(t(logZ_vec[-c(1,2),]),names=c("RB","BH","sAIS","mAIS","Boot","Adap")[-c(1,2)]); abline(h=lZ_true)
var_vec<-sapply(1:nrow(logZ_vec),function(n) var(exp(logZ_vec[n,])))
mean_vec<-sapply(1:nrow(logZ_vec),function(n) mean(exp(logZ_vec[n,])))
MSE_vec<-(mean_vec-exp(lZ_true))^2+var_vec
var_log_vec<-sapply(1:nrow(logZ_vec),function(n) var(logZ_vec[n,]))
mean_log_vec<-sapply(1:nrow(logZ_vec),function(n) mean(logZ_vec[n,]))
MSE_log_vec<-(mean_log_vec-lZ_true)^2+var_log_vec
print(round(MSE_vec/min(MSE_vec),2))
print(round(MSE_log_vec/min(MSE_log_vec),2))

rep<-100
logZ_vec<-mcmapply(function(i){
  
  KernSteps<-1;
  tmp_BH_Boot_09<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=1,beta_thresh=0.9)
  tmp_BH_SMC_09<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.9)
  tmp_BH_Boot_099<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=1,beta_thresh=0.99)
  tmp_BH_SMC_099<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.99)
  tmp_BH_Boot_0999<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=1,beta_thresh=0.999)
  tmp_BH_SMC_0999<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.999)
  
  return(c(tmp_BH_Boot_09$logZ,tmp_BH_SMC_09$logZ,tmp_BH_Boot_099$logZ,tmp_BH_SMC_099$logZ,tmp_BH_Boot_0999$logZ,tmp_BH_SMC_0999$logZ))
},1:rep)

length(beta_vec)
boxplot(t(logZ_vec),names=c("B0.9","A0.9","B0.99","A0.99","B0.999","A0.999")); abline(h=lZ_true)
var_vec<-sapply(1:nrow(logZ_vec),function(n) var(exp(logZ_vec[n,])))
mean_vec<-sapply(1:nrow(logZ_vec),function(n) mean(exp(logZ_vec[n,])))
MSE_vec<-(mean_vec-exp(lZ_true))^2+var_vec
var_log_vec<-sapply(1:nrow(logZ_vec),function(n) var(logZ_vec[n,]))
mean_log_vec<-sapply(1:nrow(logZ_vec),function(n) mean(logZ_vec[n,]))
MSE_log_vec<-(mean_log_vec-lZ_true)^2+var_log_vec
print(round(MSE_vec/min(MSE_vec),2))
print(round(MSE_log_vec/min(MSE_log_vec),2))

#save(d,lZ_true,beta_vec,logZ_vec,file="Ex4_logZ_d5.RData")
#save(d,lZ_true,logZ_vec,file="Ex4_logZ_d10_adap.RData")
load("Ex4_logZ_d10.RData")

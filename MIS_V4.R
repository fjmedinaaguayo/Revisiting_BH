library(matrixStats)
library(parallel)
library(Rcpp)
library(abind)
library(scales)
library(pbmcapply)
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
  logNum<-beta_t*logUnPi(x)+logQ_L_cpp(x,mu_ln,sigma_ln,ordrd)
  
  #tmp<-sapply(1:nrow(mu_vec_sub),function(i) beta_t*logQ_L(x,mu_vec_sub[i,]) + log(N_vec[i]) )
  #tmp<-logGamma_den_cpp(x,mu_vec_sub,N_vec,rep(sigma_prop,d),beta_t)
  logQ_vec<-logQ_L_vec_cpp(x,mu_vec_sub,sigma_vec_sub,ordrd)
  tmp<-beta_t*logQ_vec+log(N_vec)
  
  logDen<-logSumExp(tmp)
  
  out<-logNum-logDen
  
  if(is.nan(out))
    browser()
  
  return(out) 
}

logQ_L<-function(x,mu_l,sigma_l){
  
  d<-length(x)
  #out<-sum(dnorm(x,mean=mu_l,sd=sigma_prop,log=T))
  out<-logQ_L_cpp(x,mu_l,sigma_l,ordrd)
  return(out)
}

randQ_L<-function(mu_l,sigma_l){
  
  mu<-as.double(mu_l)
  
  out<-rnorm(d,mean=mu,sd=sigma_l)
  if(ordrd)
    out[1:(d/2)]<-sort(out[1:(d/2)])
  
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


MIS<-function(d,mu_vec,sigma_vec,N,plots=F,estim="all"){
  
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
                       logSumExp( sapply(1:nrow(mu_vec),function(l) logQ_L_cpp(X[,n,1],mu_vec[l,],sigma_vec[l,],ordrd)) -log(nrow(mu_vec)) ))
  
  logZ_BH<-logSumExp(logW_BH)-log(N)
  logZ_RB<-logSumExp(logW_RB)-log(N)
  
  return(list(logW_BH=logW_BH,
              logW_RB=logW_RB,
              X=X[,,1],
              L=L,
              logZ_BH=logZ_BH,
              logZ_RB=logZ_RB))
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
    logQ_vec_den<-sapply(1:N,function(n) logQ_L_vec_cpp(X[,n,t],mu_vec_sub,sigma_vec_sub,ordrd) )
    logQ_vec_num<-sapply(1:N,function(n) logQ_L_cpp(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],ordrd) )
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
                         logQ_L_cpp(X_tp1[,n],mu_vec[L_tp1[n],],sigma_vec[L_tp1[n],],ordrd) ,1:N)
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
    logQ_vec_den<-sapply(1:N,function(n) logQ_L_vec_cpp(X[,n,t],mu_vec_sub,sigma_vec_sub,ordrd) )
    logQ_vec_num<-sapply(1:N,function(n) logQ_L_cpp(X[,n,t],mu_vec[L[n,t],],sigma_vec[L[n,t],],ordrd) )
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
      
      if(log_relCESS_fn(1)>0)
        beta_adap[t+1]<-1
      else
        beta_adap[t+1]<-uniroot(log_relCESS_fn,c(beta_adap[t],1),tol=1e-100)$root
      
      if(beta_adap[t+1]==beta_adap[t])
        beta_adap[t+1] = beta_adap[t]+1e-100
      
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
physical_cores <- parallel::detectCores(logical = FALSE)
print(physical_cores)

d<-12; N<-500; KernSteps<-1; limInf<--30; limSup<-30
ordrd<-TRUE
K<-max(floor(N^(1/d)),2)
K<-max(floor((N/2)^(1/d)),2)
mu<-seq(limInf,limSup,length.out = 2*K+1)[2*(1:K)]
sigma_MCMC<-max(diff(mu)/2)/sqrt(d)
sigma_MCMC<-2.38/sqrt(d)
mu_vec<-as.matrix(expand.grid(replicate(d,mu,simplify = F)))
sigma_vec<-matrix(rep(max(diff(mu)/2),d),nrow=nrow(mu_vec),ncol=d)


n_clust<-N
K<-N
n_data<-length(y)
set.seed(123)
clust_matrix<-t(sapply(1:n_clust,function(i) sample(1:(d/2),n_data,rep=TRUE)))
emp_mus <- t(sapply(1:n_clust,function(l) sort(as.numeric(tapply(y, clust_matrix[l,], mean, na.rm = TRUE)))))
emp_sigmas <- t(sapply(1:n_clust,function(l) as.numeric(tapply(y, clust_matrix[l,], sd, na.rm = TRUE))))

if(d==2){
  emp_mus <- sapply(1:n_clust,function(l) sort(as.numeric(tapply(y, clust_matrix[l,], mean, na.rm = TRUE))))
  emp_sigmas <- sapply(1:n_clust,function(l) as.numeric(tapply(y, clust_matrix[l,], sd, na.rm = TRUE)))
  
  
}

mu_vec<-cbind(emp_mus,log(emp_sigmas))
sigma_vec<-mu_vec*0+10

source('MIS_V3_target.R')
lUnDens_list<-list(ld0,ld1,ld2,ld3)
lUnDens_list<-list(ld_mixture_model)
lWeights<-c(0)
lZ_true_vec<-c(ld0(0,"lZ",d),ld1(0,"lZ",d),ld2(0,"lZ",d),ld3(0,"lZ",d))
lZ_true<-logSumExp(lZ_true_vec+lWeights)-logSumExp(lWeights)
lZ_true <- NA

t0_start <- proc.time()
AIS_results<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.99,plots=T)
t0 <- as.numeric(proc.time() - t0_start)[3]
print(paste0("One run took: ", ceiling(t0), " secs."))
print(c(lZ_true,AIS_results$logZ))
beta_vec<-AIS_results$beta_vec


rep <- 30
all_results <- pbmclapply(1:rep, function(i){
  
  beta_thresh <- 0.99
  
  t1_start <- proc.time()
  tmp_BH_Boot <- AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,
                         alpha_thresh=1,beta_thresh=beta_thresh)
  t1 <- as.numeric(proc.time() - t1_start)[3]
  
  t2_start <- proc.time()
  tmp_BH_SMC <- AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,
                        alpha_thresh=0.5,beta_thresh=beta_thresh)
  t2 <- as.numeric(proc.time() - t2_start)[3]
  
  t3_start <- proc.time()
  tmp_BH_AIS <- AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,
                        alpha_thresh=0,beta_thresh=beta_thresh,beta_vec=tmp_BH_SMC$beta_vec)
  t3 <- as.numeric(proc.time() - t3_start)[3]
  
  t0_start <- proc.time()
  tmp_MIS<-MIS(d,mu_vec,sigma_vec,2*N*(length(tmp_BH_SMC$beta_vec)-1))
  t0 <- as.numeric(proc.time() - t0_start)[3]
  
  list(
    results = c(tmp_MIS$logZ_RB,
                tmp_MIS$logZ_BH,
                tmp_BH_AIS$logZ, 
                tmp_BH_Boot$logZ, 
                tmp_BH_SMC$logZ),
    times = c(MIS = t0,
              BH_AIS = t3, 
              BH_Boot = t1, 
              BH_SMC = t2),
    n_dists = c(length(tmp_BH_AIS$beta_vec),
                length(tmp_BH_Boot$beta_vec),
                length(tmp_BH_SMC$beta_vec))
  )
}, mc.cores = max(physical_cores-1,1))

# Extract results and times
logZ_vec <- sapply(all_results, function(x) x$results)
time_vec <- sapply(all_results, function(x) x$times)
n_dists_vec <- sapply(all_results, function(x) x$n_dists)

boxplot(t(logZ_vec),names=c("RB","BH","mAIS","Boot","Adap")); abline(h=lZ_true)
rowMeans(time_vec)
boxplot(t(time_vec),names=c("MIS","mAIS","Boot","Adap"));
boxplot(t(n_dists_vec),names=c("mAIS","Boot","Adap"));

save(logZ_vec,time_vec,n_dists_vec,rep,d,N,K,file="./TEMP_mixmodel_b99_N500_d12_v2.RData")
# load(file="./mixmodel_N500_d14.RData")

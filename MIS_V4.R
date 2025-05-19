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
                       logSumExp( sapply(1:nrow(mu_vec),function(l) logQ_L_cpp(X[,n,1],mu_vec[l,],sigma_vec[l,],ordrd)) -log(nrow(mu_vec)) ))
  
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
  
  X<-array(0,dim=c(d,N,1))
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
    X<-abind(X, array(0, replace(dim(X), 3, 1)), along = 3)
    
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


d<-3; N<-1000; KernSteps<-1; limInf<--10; limSup<-10
K<-max(floor(N^(1/d)),2)
K<-max(floor((10*N)^(1/d)),2)
ordrd<-FALSE
mu<-seq(limInf,limSup,length.out = 2*K+1)[2*(1:K)]
sigma_MCMC<-max(diff(mu)/2)/sqrt(d)
sigma_MCMC<-2.38/sqrt(d)
mu_vec<-as.matrix(expand.grid(replicate(d,mu,simplify = F)))
sigma_vec<-matrix(rep(max(diff(mu)/2),d),nrow=nrow(mu_vec),ncol=d)

source('./MIS_V3_target.R')
lUnDens_list<-list(ld0,ld1,ld2,ld3)
lUnDens_list<-list(ld_mixture_model)
lUnDens_list<-list(ld_my_rosenbrock)
lWeights<-c(0)
lZ_true_vec<-c(ld0(0,"lZ",d),ld1(0,"lZ",d),ld2(0,"lZ",d),ld3(0,"lZ",d))
lZ_true<-logSumExp(lZ_true_vec+lWeights)-logSumExp(lWeights)
lZ_true <- ld_my_rosenbrock(rep(0,d),"lZ",d)

AIS_results<-AIS_MIS(d,mu_vec,sigma_vec,N,KernSteps,sigma_MCMC,alpha_thresh=0.5,beta_thresh=0.9,plots=T)
print(c(lZ_true,AIS_results$logZ))
beta_vec<-AIS_results$beta_vec

i1<-1; i2<-2
x1<-seq(limInf,limSup,.1); x2<-x1
targetUnDens<-function(x,y){
  
  theta<-rep(0,d)
  theta[c(i1,i2)]<-c(x,y)
  
  res<-logTargetUnDens(theta,lUnDens_list,lWeights)
  
  return(res)
}
z<-outer(x1,x2,FUN=Vectorize(targetUnDens))

plot(matrix(1:ncol(AIS_results$L),nrow=N,ncol=ncol(AIS_results$L),byrow = T),AIS_results$L)

library(ggplot2)
library(tidyr)
library(dplyr)

# Prepare the data
plot_data <- data.frame(
  x = as.vector(AIS_results$X[1,1:1000,]),  # First component
  y = as.vector(AIS_results$X[2,1:1000,]),  # Second component
  z = as.vector(AIS_results$X[3,1:1000,]),  # Second component
  time = rep(1:dim(AIS_results$X)[3], each = 1000)
)


# Create the ggplot
plot_data %>% 
  filter(time %in% c(1,5,10,15,20,30,ncol(AIS_results$L))) %>%
  ggplot(aes(x = x, y = y, color = factor(time))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    x = "x[1]",
    y = "x[2]",
    color = "s"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

plot_data %>% 
  filter(time %in% c(1, 5, 10, 15, 20, 30)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~ time, scales = "free", 
             labeller = labeller(time = function(t) paste0("s=", t))) +
  labs(
    x = expression(x[1]),
    y = expression(x[2])
  ) +
  theme_minimal()

ggsave("distn_particles_ex2_bw.png", width = 10, height = 8, dpi = 300, units = "in")

# Create the ggplot
plot_data %>% 
  filter(time %in% c(1,5,10,15,20,30,ncol(AIS_results$L))) %>%
  ggplot(aes(x = x, y = z, color = factor(time))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    x = "x[2]",
    y = "x[3]",
    color = "s"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Create the ggplot
plot_data %>% 
  filter(time %in% c(1,5,10,15,20,30,ncol(AIS_results$L))) %>%
  ggplot(aes(x = y, y = z, color = factor(time))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    x = "x[2]",
    y = "x[3]",
    color = "s"
  ) +
  theme_minimal() +
  theme(legend.position = "right")



library(ggridges)
library(dplyr)

# Prepare the data
# Assuming you want to plot the first component across time
plot_data <- data.frame(
  value = as.vector(AIS_results$X[1,1:1000,]),  # First component
  time = factor(rep(1:dim(AIS_results$X)[3], each = 1000), 
                levels = 1:dim(AIS_results$X)[3])  # Ensure correct ordering
)

# Create the joyplot
plot_data %>% 
  filter(time %in% c(1,2,4,8,16,32)) %>%
  ggplot(aes(x = value, y = time, fill = time)) +
  geom_density_ridges(
    fill = "black",
    scale = 0.9,  # Adjust overlap of ridges
    color = "black",  # Outline color
  ) +
  labs(
    x = "x",
    y = "s"
  ) +
  theme_ridges() +
  theme(legend.position = "none")

# Create the joyplot in the style of Unknown Pleasures album cover
plot_data %>% 
  ggplot(aes(x = value, y = as.factor(time), height = ..density..)) +
  geom_density_ridges_gradient(
    fill = "black",
    color = "white",
    scale = 10
  ) +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

rep <- 30
all_results <- pbmclapply(1:rep, function(i){
  
  beta_thresh <- 0.9
  
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
  #tmp_MIS<-MIS(d,mu_vec,sigma_vec,N)
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
}, mc.cores = physical_cores-1)

# Extract results and times
logZ_vec <- sapply(all_results, function(x) x$results)
time_vec <- sapply(all_results, function(x) x$times)
n_dists_vec <- sapply(all_results, function(x) x$n_dists)

boxplot(t(logZ_vec),names=c("RB","BH","AIS","Boot","SMC")); abline(h=lZ_true)
boxplot(t(logZ_vec)[,1:2],names=c("RB","BH")); abline(h=lZ_true)
boxplot(t(logZ_vec)[,4:5],names=c("Boot","SMC")); abline(h=lZ_true)
rowMeans(time_vec)
boxplot(t(time_vec),names=c("MIS","AIS","Boot","Adap"));
boxplot(t(n_dists_vec),names=c("mAIS","Boot","Adap"));

# save(all_results,file="./rosenbrock_3d_4paper.RData")

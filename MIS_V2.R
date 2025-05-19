library(matrixStats)
library(mvtnorm)
library(parallel)

logTargetUnDens<-function(theta,lUnDens_list,lWeights){
  
  d<-length(theta)
  nDens<-length(lUnDens_list)
  lNormWeights<-lWeights-logSumExp(lWeights)
  
  lDens_eval<-sapply(1:nDens, function(i) lUnDens_list[[i]](theta,d=d))
  
  val<-logSumExp(lDens_eval+lNormWeights)
  
  return(val)
}

randTarget<-function(N,lUnDens_list,lWeights){
  
  nDens<-length(lUnDens_list)
  lNormWeights<-lWeights-logSumExp(lWeights)
  
  idDist<-sample.int(nDens,N,rep=T,prob=exp(lNormWeights))
  
  draws<-t(sapply(idDist, function(i) lUnDens_list[[i]](theta,type="draw")))
  
  return(draws)
}

logProposalDens<-function(theta,type,K=1,mu_vec=numeric(),mu_list=list(),Sigma_list=list(),N_vec=numeric()){
  
  d<-length(theta)
  
  K_mu<-length(mu_list)
  K_Sigma<-length(Sigma_list)
  
  if(K_mu>0 & K_Sigma>0 & K_mu==K_Sigma){
    
    lDens_vec<-sapply(1:K_mu,function(i) dmvnorm(theta,mean=as.double(mu_list[[i]]),sigma=Sigma_list[[i]],log=T))
    
    if(length(N_vec)==K_mu)
      lDens_vec<-lDens_vec+log(N_vec)
    
    return(logSumExp(lDens_vec))
  }
  
  sigma_prop<-limSupp/(q*K)
  Cov_prop<-diag(rep(sigma_prop^2,d))
  
  if(type=="IS_Gauss" & K==1){
    
    return(sum(dnorm(theta,rep(0,d),rep(sigma_prop,d),log=T)))
    #return(dmvnorm(theta,sigma=Cov_prop,log=T))
  }
  else if(type=="IS_unif" & K==1){
    
    return(-d*log(2*limSupp))
  }
  else if(type=="BH" & length(mu_vec)>0){
    
    lDens_vec<-sapply(1:nrow(mu_vec), function(i) sum(dnorm(theta,as.double(mu_vec[i,]),rep(sigma_prop,d),log=T)))
    #lDens_vec<-sapply(labels,function(i) dmvnorm(theta,mean=as.double(mu_vec[i,]),sigma=Cov_prop,log=T))
    
    return(logSumExp(lDens_vec))
  }
  else if(type=="marg" & length(mu_vec)>0){
    
    lDens_vec<-sapply(1:nrow(mu_vec), function(i) sum(dnorm(theta,as.double(mu_vec[i,]),rep(sigma_prop,d),log=T)))
    #lDens_vec<-sapply(1:nrow(mu_vec), function(i) dmvnorm(theta,mean=as.double(mu_vec[i,]),sigma=Cov_prop,log=T))
    
    return(logSumExp(lDens_vec)-d*log(K))
  }
}

randProposal<-function(N,type,K=1,d=2,mu_list=list(),Sigma_list=list()){
  
  K_mu<-length(mu_list)
  K_Sigma<-length(Sigma_list)
  
  if(K_mu>0 & K_Sigma>0 & K_mu==K_Sigma){
    
    label<-sample.int(K_mu,N,replace=T)
    theta<-t(sapply(label,function(i) rmvnorm(1,as.double(mu_list[[i]]),sigma=Sigma_list[[i]] )))
    
    return(list(label,theta))
  }
  
  mu_prop_min<--limSupp+limSupp/K
  mu_prop_max<-limSupp-limSupp/K
  mu<-seq(mu_prop_min,mu_prop_max,length.out=K)
  mu_vec<-expand.grid(replicate(d,mu,simplify = F))
  sigma_prop<-limSupp/(q*K)
  Cov_prop<-diag(rep(sigma_prop^2,d))
  
  if(type=="IS_Gauss" & K==1){
    
    label<-rep(1,N)
    theta<-rmvnorm(N,rep(0,d),sigma=Cov_prop)
    
    return(list(label,theta))
  }
  if(type=="IS_unif" & K==1){
    
    label<-rep(1,N)
    theta<-matrix(runif(N*d,-limSupp,limSupp),ncol=d)
    
    return(list(label,theta))
  }
  else if(type=="BH" | type=="marg"){
    
    label<-sample.int(nrow(mu_vec),N,replace = T)
    theta<-t(sapply(label,function(i) rmvnorm(1,as.double(mu_vec[i,]),sigma=Cov_prop)))
    
    return(list(label,theta))
  }
}

ld0<-function(theta,type="lUnDen",d){

  mu<-rep(0,d)
  diagVec<-rep(10,d)
  
  Sigma<-diag(diagVec)
  Sigma_inv<-diag(1/diagVec)
  
  if(type=="lZ")
    return(d/2*log(2*pi)+0.5*log(abs(det(Sigma))))
  if(type=="draw")
    return(rmvnorm(1,mu,Sigma))
  else if(type=="lUnDen")
    return(-0.5*sum((theta-mu)%*%Sigma_inv%*%(theta-mu)))
}

ld1<-function(theta,type="lUnDen",d){
  
  Transf<-function(x,b){
    
    x[2]<-x[2]+b*x[1]^2-100*b
    
    return(x)
  }
  
  mu<-rep(0,d)
  diagVec<-rep(1,d)
  diagVec[1]<-diagVec[1]**100
  Sigma<-diag(diagVec)
  Sigma_inv<-diag(1/diagVec)
  
  if(type=="lZ")
    return(d/2*log(2*pi)+0.5*log(abs(det(Sigma))))
  
  b<-0.03
  theta<-Transf(theta,b)
  
  if(type=="draw")
    return(NA)
  else if(type=="lUnDen")
    return(-0.5*sum((theta-mu)%*%Sigma_inv%*%(theta-mu)))
}

ld2<-function(theta,type="lUnDen",d){
  
  mu<-rep(-10,d)
  Sigma<-diag(rep(1,d))
  Sigma[1,2]<-2
  Sigma[2,1]<-2
  Sigma[2,2]<-5
  Sigma_inv<-solve(Sigma)
  
  if(type=="draw")
    return(rmvnorm(1,mu,Sigma))
  else if(type=="lZ")
    return(d/2*log(2*pi)+0.5*log(abs(det(Sigma))))
  else if(type=="lUnDen")
    return(-0.5*sum((theta-mu)%*%Sigma_inv%*%(theta-mu)))
}

ld3<-function(theta,type="lUnDen",d){
  
  nu<-3
  mu<-rep(5,d)
  Sigma<-diag(rep(1,d))
  Sigma_inv<-solve(Sigma)
  
  if(type=="draw")
    return(rmvt(1,sigma=Sigma,df=nu,delta=mu))
  else if(type=="lZ")
    return(lgamma(nu/2)-lgamma(0.5*(nu+d))+d/2*log(nu*pi)+0.5*log(abs(det(Sigma))))
  else if(type=="lUnDen")
    return(-0.5*(nu+d)*log(1+1/nu*sum((theta-mu)%*%Sigma_inv%*%(theta-mu))))
}

lUnDens_list<-list(ld1)
lWeights<-c(0)
lZ_true_vec<-c(ld1(0,"lZ",d=d))

# lUnDens_list<-list(ld1,ld2)
# lWeights<-c(0,0)
# lZ_true_vec<-c(ld1(0,"lZ"),ld2(0,"lZ"))

# lUnDens_list<-list(ld1)
# lWeights<-c(0)
# lZ_true_vec<-c(ld1(0,T))

lZ_true<-logSumExp(lZ_true_vec+lWeights)-logSumExp(lWeights)

theta<-c(0,0,0)

logTargetUnDens(theta,lUnDens_list,lWeights)
logProposalDens(theta,type='IS_Gauss',)
logProposalDens(theta,type='IS_unif')
logProposalDens(theta,type='BH',K=10)
logProposalDens(theta,type='marg',K=10)


logZ_IS<-function(N,rep,type){
  
  lZ_rep<-mcmapply(function(j){
    
    sample<-randProposal(N,type,d=d)
    theta<-sample[[2]]
    
    lW<-sapply(1:N,function(i){
      logTargetUnDens(theta[i,],lUnDens_list,lWeights)-logProposalDens(theta[i,],type)
      })
    lZ<-logSumExp(lW)-log(N)
    
    return(lZ)
  },1:rep)
  
  return(lZ_rep)
}

logZ_BH<-function(N,rep,K){
  
  mu_prop_min<--limSupp+limSupp/K
  mu_prop_max<-limSupp-limSupp/K
  mu<-seq(mu_prop_min,mu_prop_max,length.out=K)
  mu_vec<-expand.grid(replicate(d,mu,simplify = F))
    
  lZ_rep<-mcmapply(function(j){
    
    sample<-randProposal(N,"BH",K,d)
    labls<-sample[[1]]
    theta<-sample[[2]]
    
    lW<-sapply(1:N,function(i){
      logTargetUnDens(theta[i,],lUnDens_list,lWeights)-logProposalDens(theta[i,],"BH",K,mu_vec[labls,])
    })
    lZ<-logSumExp(lW)
    
    return(lZ)
  },1:rep)
  
  return(lZ_rep)
}

logZ_BH2<-function(N,rep,mu_list,Sigma_list){
  
  lZ_rep<-mcmapply(function(j){
    
    sample<-randProposal(N,"BH",mu_list=mu_list,Sigma_list=Sigma_list)
    labls<-sample[[1]]
    theta<-sample[[2]]
    
    N_vec<-table(labls)
    ind_N_vec<-as.double(names(N_vec))
    
    lW<-sapply(1:N,function(i){
      logTargetUnDens(theta[i,],lUnDens_list,lWeights)-
        logProposalDens(theta[i,],"BH",mu_list=mu_list[ind_N_vec],Sigma_list=Sigma_list[ind_N_vec],N_vec=N_vec)
    })
    lZ<-logSumExp(lW)
    
    return(lZ)
  },1:rep)
  
  return(lZ_rep)
}

logZ_marg2<-function(N,rep,mu_list,Sigma_list){
  
  K<-length(mu_list)
  
  lZ_rep<-mcmapply(function(j){
    
    sample<-randProposal(N,"marg",mu_list=mu_list,Sigma_list=Sigma_list)
    labls<-sample[[1]]
    theta<-sample[[2]]
    
    lW<-sapply(1:N,function(i){
      logTargetUnDens(theta[i,],lUnDens_list,lWeights)-
        logProposalDens(theta[i,],"marg",mu_list=mu_list,Sigma_list=Sigma_list)+log(K)
    })
    lZ<-logSumExp(lW)-log(N)
    
    return(lZ)
  },1:rep)
  
  return(lZ_rep)
}

logZ_marg<-function(N,rep,K){
  
  mu_prop_min<--limSupp+limSupp/K
  mu_prop_max<-limSupp-limSupp/K
  mu<-seq(mu_prop_min,mu_prop_max,length.out=K)
  mu_vec<-expand.grid(replicate(d,mu,simplify = F))
  
  lZ_rep<-sapply(1:rep,function(j){
    
    sample<-randProposal(N,"marg",K,d)
    labels<-sample[[1]]
    theta<-sample[[2]]
    
    lW<-sapply(1:N,function(i){
      logTargetUnDens(theta[i,],lUnDens_list,lWeights)-logProposalDens(theta[i,],"marg",K,mu_vec)
    })
    lZ<-logSumExp(lW)-log(N)
    
    return(lZ)
  })
  
  return(lZ_rep)
}

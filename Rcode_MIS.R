library(matrixStats)
library(rmutil) #for betabinomial
library(parallel)
library(ggplot2)
library(reshape2)

log_un_eta<-function(index,theta_vec,label_vec,gamma,scheme){
  
  if(scheme=="PG"){
    
    log_den<-gamma*(log_un_pi(theta_vec[index])-logSumExp(log_q_cond(theta_vec[index],label_vec)))+
      sum(log_q(theta_vec,label_vec))-(1-gamma)*log(length(theta_vec))
  }
  else if(scheme=="SG"){
    
    log_den<-gamma*log_un_pi(theta_vec[index])-logSumExp(gamma*log_q_cond(theta_vec[index],label_vec))+
      sum(log_q(theta_vec,label_vec))
  }
  
  return(log_den)
}

log_un_pi<-function(theta){
  
  return(dt(theta,df_pi,mu_pi,log=T)+lZ_true)
}

log_alpha<-function(l){
  
  if(s==Inf)
    log_den_l<-dbinom(l-1,K-1,m,log=T)
  else if(m==0.5 & s==2)
    log_den_l<--log(K)
  else
    log_den_l<-lbeta(l-1+s*m,K-1-(l-1)+s*(1-m)) + lchoose(K-1,l-1) -lbeta(s*m,s*(1-m))
    #log_den_l<-dbetabinom(l-1,K-1,m,s,log=T)
    
    
  return(log_den_l)
}

log_q_cond<-function(theta,l){
  
  mu_prop<-(mu_prop_max-mu_prop_min)*(l-1)/(K-1)+mu_prop_min
  log_den_theta<-dnorm(theta,mu_prop,sigma_prop,log=T)
  
  return(log_den_theta)
}

log_q<-function(theta,l){
  
  return(log_q_cond(theta,l)+log_alpha(l))
}

log_q_marg<-function(theta){
  
  log_den<-logSumExp(sapply(1:K, function(k) log_q(theta,k)))
  
  return(log_den)
}

r_q<-function(N){ #Sample is a Matrix of N rows, 1st column is label
  
  if(s==Inf)
    l<-rbinom(N,K-1,rep(m,N))+1
  else
    l<-rbinom(N,K-1,rbeta(N,s*m,s*(1-m)))+1
  
  mu_prop<-(mu_prop_max-mu_prop_min)*(l-1)/(K-1)+mu_prop_min
  theta<-rnorm(N,mu_prop,sigma_prop)
  
  return(matrix(c(l,theta),nrow=N))
}

r_q_cond<-function(l){
  
  mu_prop<-(mu_prop_max-mu_prop_min)*(l-1)/(K-1)+mu_prop_min
  theta<-rnorm(length(l),mu_prop,sigma_prop)
  
  return(theta)
}

altRep<-function(Particles,lW=NA){ #Sample is a List of size K_eff, each element also a list of size 2: 1st element is label, 2nd is vector of values
  
  labels<-as.numeric(names(table(Particles[,1])))
  if(max(is.na(lW))){
    
    result<-lapply(labels, function(l) list(l,Particles[which(Particles[,1]==l),2]))
  }
  else{
    
    result<-lapply(labels, function(l){
      
      ind_vec<-which(Particles[,1]==l)
      PartsNlW<-cbind(Particles,lW[])[ind_vec,2:3]
      
      return(list(c(l,PartsNlW)))
    })
  }
  
  return(result)
}

obtain_Nj<-function(Particles_alt){ #Compute vector of number of values for each label
  
  K_eff<-length(Particles_alt)
  N_j<-sapply(1:K_eff, function(k) length(Particles_alt[[k]][[2]]) )
  
  return(N_j)
}

IS_opt_AIS<-function(gamma_vec,Particles,MCMCtype="none",iters=1,verbose=F){
  
  N<-nrow(Particles)
  AIS_steps<-length(gamma_vec)
  
  lW<-log_un_pi(Particles[,2])-sapply(1:N, function(i) log_q_marg(Particles[i,2]))
  
  lZ<-logSumExp(lW)-log(N)
  
  return(lZ)
}

BH_AIS_std<-function(gamma_vec,Particles,Particles_alt,scheme,MCMCtype="none",iters=1,verbose=F,sigma_prop_MCMC){
  
  N<-nrow(Particles)
  K_eff<-length(Particles_alt)
  N_j<-obtain_Nj(Particles_alt)
  AIS_steps<-length(gamma_vec)

  lW<-logSumExp(sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[1],scheme))) -
    logSumExp(sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],0,scheme)))
  
  if(AIS_steps==1)
    return(lW)
  
  for(t in 1:(AIS_steps-1)){
    
    # indices<-sample(1:N,N,rep=T,prob=exp(sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[1],scheme))))
    # print(cbind(Particles,Particles[indices,]))
    # Particles<-Particles[indices,]
    
    if(MCMCtype=="Gibbs")
      Particles<-Gibbsmove_BH(Particles,iters,gamma_vec[t],scheme,verbose)
    else if(MCMCtype=="MCMC")
      Particles<-MCMCmoveAll_BH(Particles,iters,gamma_vec[t],scheme,verbose,sigma_prop_MCMC)
    
    lW[t+1]<-logSumExp(sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[t+1],scheme))) -
      logSumExp(sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[t],scheme)))

  }
  
  lZ<-sum(lW)
  
  return(lZ)
}

BH_AIS_mod<-function(gamma_vec,Particles,Particles_alt,scheme,MCMCtype="none",iters=1,verbose=F,sigma_prop_MCMC){
  
  N<-nrow(Particles)
  K_eff<-length(Particles_alt)
  N_j<-obtain_Nj(Particles_alt)
  AIS_steps<-length(gamma_vec)
  
  lW<-sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[1],scheme)) -
    logSumExp(sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],0,scheme)))
  
  if(AIS_steps==1)
    return(logSumExp(lW))
  
  for(t in 1:(AIS_steps-1)){
    
    if(MCMCtype=="MCMC"){
      
      Particles[,2]<-sapply(1:N, function(n) MCMCmove_BH(n, Particles,iters,gamma_vec[t],scheme,verbose,sigma_prop_MCMC))
    }
    
    if(verbose){
      
      #lines(density(Particles[,2]),col="red",lwd=1,lty=3)
    }
    
    lW<-lW+sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[t+1],scheme)) -
      sapply(1:N, function(i) log_un_eta(i,Particles[,2],Particles[,1],gamma_vec[t],scheme))
  }
  
  lZ<-logSumExp(lW)
  
  return(lZ)
}

Gibbsmove_BH<-function(Particles,iters,gamma,scheme,verbose){
  
  N<-nrow(Particles)
  log_probs<-matrix(0,nrow=N,ncol=iters)
  temp<-Particles[,2]
  temp2<-numeric(0)
  
  for(i in 1:iters){
    
    log_probs[,i]<-sapply(1:N, function(k){
      
      log_un_eta(k,Particles[,2],Particles[,1],gamma,scheme)
    })
    
    log_probs[,i]<-log_probs[,i]-max(log_probs[,i])
    n<-sample(1:N,1,prob=exp(log_probs[,i]))
    temp2[i]<-n
    
    Particles[-n,2]<-r_q_cond(Particles[-n,1])
    temp<-cbind(temp,Particles[,2])
  }
  
  if(verbose){
    
    #boxplot(exp(log_probs))
    #lines(density(Particles[,2]),col="red",lwd=1,lty=3)
    # m<-ggplot(as.data.frame(temp2), aes(x=temp2)) + 
    #   geom_histogram(color="black", fill="white",binwidth = 1)
    # print(m)
    # df<-as.data.frame(t(temp))
    # m <- ggplot(df, aes(x = V1, y = V2)) +
    #   geom_point() + geom_density_2d()
    # print(m)
    # m <- ggplot(melt(df), aes(value, col=variable)) +
    #   geom_density()
    # print(m)
  }
  
  return(Particles)
}

MCMCmove_BH<-function(n,Particles,iters,gamma,scheme,verbose,sigma_prop_MCMC,GF=0,log_rho_vec=NA){
  
  N<-nrow(Particles)
  accept<-0
  
  for(i in 1:iters){
    
    Particles_star<-Particles
    Particles_star[n,2]<-Particles_star[n,2]+rnorm(1,0,sigma_prop_MCMC)
    
    if(GF==0){
      lratio<-log_un_eta(n,Particles_star[,2],Particles_star[,1],gamma,scheme)-
        (log_un_eta(n,Particles[,2],Particles[,1],gamma,scheme))
    }
    else{
      lratio<-log_un_eta_GF(n,Particles_star[,2],Particles_star[,1],gamma,scheme,log_rho_vec)-
        (log_un_eta_GF(n,Particles[,2],Particles[,1],gamma,scheme,log_rho_vec))
    }
    
    if(log(runif(1))<lratio){
      
      accept<-accept+1
      n0<-n
      Particles<-Particles_star
    }
  }
  
  if(verbose){
    
    print(paste0("Acceptance rate: ", accept/iters))
  }
  
  return(Particles[n,2])
}

MCMCmoveAll_BH<-function(Particles,iters,gamma,scheme,verbose,sigma_prop_MCMC){
  
  N<-nrow(Particles)
  accept<-0
  
  for(i in 1:iters){
    
    Particles_star<-Particles
    Particles_star[,2]<-Particles_star[,2]+rnorm(N,0,sigma_prop_MCMC/sqrt(N))
    
    lratio<-logSumExp(sapply(1:N, function(n) log_un_eta(n,Particles_star[,2],Particles_star[,1],gamma,scheme)))-
      logSumExp(sapply(1:N, function(n) log_un_eta(n,Particles[,2],Particles[,1],gamma,scheme)))
    
    if(log(runif(1))<lratio){
      
      accept<-accept+1
      Particles<-Particles_star
    }
  }
  
  if(verbose){
    
    print(paste0("Acceptance rate: ", accept/iters))
  }
  
  return(Particles)
}


#### Intractable proposal

IS_extended<-function(Particles){
  
  N<-nrow(Particles)
  
  lZ<-logSumExp(log_un_pi(Particles[,2])-log(K)-log_q(Particles[,2],Particles[,1]))-log(N)
  
  return(lZ)
}


IS_opt<-function(Particles){
  
  N<-nrow(Particles)
  
  lZ<-logSumExp(log_un_pi(Particles[,2])-sapply(1:N, function(i) logSumExp(log_q(Particles[i,2],1:K))))-log(N)
  
  return(lZ)
}

compute_lZ_G<-function(Particles,Particles_alt){ #Optimal combination of estimators using direct estimate of Sigma^-1
  
  N<-nrow(Particles)
  lZ<-IS_extended(Particles)
  K_eff<-length(Particles_alt)
  
  temp<-sapply(1:K_eff, function(j){
    
    N_j<-length(Particles_alt[[j]][[2]])
    lW_j<-sapply(1:N_j, function(i) log_un_pi(Particles_alt[[j]][[2]][i])-
                   log_q(Particles_alt[[j]][[2]][i],Particles_alt[[j]][[1]]))
    lZ_j<-logSumExp(lW_j-log(N))
    
    log_num<-sapply(1:N, function(i) 2*log_un_pi(Particles[i,2])-log(K)-
                      log_q(Particles[i,2],Particles[i,1])-log_q(Particles[i,2],Particles_alt[[j]][[1]]))
    log_diag_j<-logSumExp(log_num)-log(N)-2*lZ
    
    return(c(lZ_j,log_diag_j,N_j))
  })
  
  lZ_j<-temp[1,]
  log_diag<-temp[2,]
  N_j<-temp[3,]
  
  e<-rep(1,K_eff)
  l_min<-min(log_diag)
  
  A_0_inv<-diag(exp(-(log_diag-l_min))); B_0<-t(e)%*%A_0_inv
  Sigma_0_inv_dir<-A_0_inv*(1+exp(logSumExp(-log_diag)))-exp(-l_min)*t(B_0)%*%B_0
  lv_j_dir_0_unnorm<-log(t(e)%*%Sigma_0_inv_dir)
  lv_j_dir_0<-lv_j_dir_0_unnorm-logSumExp(lv_j_dir_0_unnorm)
  
  lZ_G<-logSumExp(lv_j_dir_0+lZ_j)
  
  return(lZ_G)
}

compute_lZ_A1<-function(Particles,Particles_alt){ #Approximation 1 to random BH
  
  N<-nrow(Particles)
  K_eff<-length(Particles_alt)
  N_j<-obtain_Nj(Particles_alt)
  
  lW_A1<-sapply(1:N, function(i){
    
    ind_Nj<-which(sapply(1:K_eff, function(j){
      
      return((Particles[i,1]==Particles_alt[[j]][[1]])*1)
    })==1)
    
    return(log_un_pi(Particles[i,2])-log_q(Particles[i,2],Particles[i,1])+log(K^-1+N_j[ind_Nj]-1)-log(N))
  })
  
  lZ_A1<-logSumExp(lW_A1)-log(N)
  
  w<-exp(lW_A1-logSumExp(lW_A1))
  R1<-sum(w*Particles[,2])
  R2<-sum(w*Particles[,2]^2)
  C2<-sum(w*(Particles[,2]-R1)^2)
  
  return(list(lZ_A1,R1,R2,C2))
}

compute_lZ_A2<-function(Particles,Particles_alt){ #Approximation 2 to random BH
  
  N<-nrow(Particles)
  K_eff<-length(Particles_alt)
  N_j<-obtain_Nj(Particles_alt)
  
  lW_A2<-sapply(1:N, function(i){
    
    log_den<-sapply(1:K_eff, function(j){
      log_q(Particles[i,2],Particles_alt[[j]][[1]]) + log(N_j[j])
    })
    
    ind_Nj<-which(sapply(1:K_eff, function(j){
      
      return((Particles[i,1]==Particles_alt[[j]][[1]])*1)
    })==1)
    
    return(log_un_pi(Particles[i,2])-logSumExp(log_den)+log(K^-1+N_j[ind_Nj]-1)-log(N))
  })
  
  lZ_A2<-logSumExp(lW_A2)
  
  w<-exp(lW_A2-logSumExp(lW_A2))
  R1<-sum(w*Particles[,2])
  R2<-sum(w*Particles[,2]^2)
  C2<-sum(w*(Particles[,2]-R1)^2)
  
  return(list(lZ_A2,R1,R2,C2))
}

log_psi<-function(n,theta,label_vec){
  
  log_q(theta,label_vec[n])
}

log_rho_all<-function(Particles,Particles_alt){
  
  N<-nrow(Particles)
  K_eff<-length(Particles_alt)
  N_j<-obtain_Nj(Particles_alt)
  
  log_rho<-sapply(1:N, function(i){
    
    ind_Nj<-which(sapply(1:K_eff, function(j){
      
      return((Particles[i,1]==Particles_alt[[j]][[1]])*1)
    })==1)
    
    return(log(K^-1+N_j[ind_Nj]-1)-log(N))
  })
  
  return(log_rho)
}

log_un_eta_GF<-function(index,theta_vec,label_vec,gamma,scheme,log_rho_vec){
  
  N<-length(theta_vec)
  if(scheme=="PG"){
    
    log_den<-gamma*(log_un_pi(theta_vec[index])+
                      log_psi(index,theta_vec[index],label_vec)+
                      log_rho_vec[index]-
                      logSumExp(log_psi(1:N,theta_vec[index],label_vec))-
                      log_q(theta_vec[index],label_vec[index]))+
      sum(log_q(theta_vec,label_vec))-(1-gamma)*log(N)
  }
  else if(scheme=="SG"){
    
    log_den<-gamma*(log_un_pi(theta_vec[index])+
                      log_psi(index,theta_vec[index],label_vec)+
                      log_rho_vec[index]-
                      log_q(theta_vec[index],label_vec[index]))-
      logSumExp(gamma*log_psi(1:N,theta_vec[index],label_vec))+
      sum(log_q(theta_vec,label_vec))
  }
  
  return(log_den)
}

GF_AIS_mod<-function(gamma_vec,Particles,Particles_alt,scheme,MCMCtype="none",iters=1,verbose=F,sigma_prop_MCMC){
  
  N<-nrow(Particles)
  K_eff<-length(Particles_alt)
  N_j<-obtain_Nj(Particles_alt)
  AIS_steps<-length(gamma_vec)
  log_rho_vec<-log_rho_all(Particles,Particles_alt)
  
  lW<-sapply(1:N, function(i) log_un_eta_GF(i,Particles[,2],Particles[,1],gamma_vec[1],scheme,log_rho_vec)) -
    logSumExp(sapply(1:N, function(i) log_un_eta_GF(i,Particles[,2],Particles[,1],0,scheme,log_rho_vec)))
  
  if(AIS_steps==1)
    return(logSumExp(lW))
  
  for(t in 1:(AIS_steps-1)){
    
    if(MCMCtype=="MCMC"){
      
      Particles[,2]<-sapply(1:N, function(n) MCMCmove_BH(n, Particles,iters,gamma_vec[t],scheme,verbose,sigma_prop_MCMC,GF=1,log_rho_vec=log_rho_vec))
    }
    
    if(verbose){
      
      #lines(density(Particles[,2]),col="red",lwd=1,lty=3)
    }
    
    lW<-lW+sapply(1:N, function(i) log_un_eta_GF(i,Particles[,2],Particles[,1],gamma_vec[t+1],scheme,log_rho_vec)) -
      sapply(1:N, function(i) log_un_eta_GF(i,Particles[,2],Particles[,1],gamma_vec[t],scheme,log_rho_vec))
  }
  
  lZ<-logSumExp(lW)
  
  return(lZ)
}

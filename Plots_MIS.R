#Target parameters
mu_pi<-0; df_pi<-Inf; lZ_true<-log(1)

#Proposals parameters
K<-3e+0 #Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-1
m<-0.5; s<-2 #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0
theta<-seq(mu_prop_min-5*sigma_prop,mu_prop_max+5*sigma_prop,.1)
plot_path<-"~/Dropbox/Reading/felipe/MIS/Plots/"

#Other parameters
N<-500
gamma_vec<-(seq(0,1,l=2)^(1))[-1]
rep<-50
print(detectCores())
options(mc.cores = detectCores()-1)
getOption("mc.cores",2L)
listEstims<-list()

# Plots for figures 1 and 2
for(K in c(3,300,30000)){
  for(m in c(0.20,0.35,0.5)){
    for(s in c(2,20,Inf)){
      
      Particles<-r_q(N)
      Particles_alt<-altRep(Particles)
      
      df_temp<-data.frame(chain=rep("Proposal",N),
                          iteration=(1:N),
                          x=Particles[,2])
      
      p1<-ggplot(df_temp, aes(x,color=chain, fill=chain))+
        stat_function(data=data.frame(x = c(-25,25), chain=rep("Target",2)), mapping=aes(x), fun = function(x) exp(log_un_pi(x)-lZ_true), geom="area")+
        geom_area(stat="density",position="identity")+
        geom_line(stat="density", size=1)+
        scale_fill_manual(values=c(NA,"lightgray"))+
        scale_color_manual(values=c("black",NA))+
        coord_cartesian(xlim=c(-25,25), ylim=c(0,0.5))+
        ylab("density")+
        theme_classic()+
        theme(legend.title=element_blank(),legend.position="none")
      #setEPS()
      #postscript(file=paste0(plot_path,"density_K",K,"_m",m,"_s",s,".eps"),width=4, height=4)
      print(p1)
      dev.copy(pdf,paste0(plot_path,"density_K",K,"_m",m,"_s",s,".pdf"),width=2, height=1.7)
      dev.off()
      
      lZ_rep<-mcmapply(function(r){
        
        Particles<-r_q(N)
        Particles_alt<-altRep(Particles)
        K_eff<-length(Particles_alt)
        lZ_BH<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG")
        Particles<-r_q(ceiling(N*K_eff/K))
        lZ_Opt<-IS_opt_AIS(gamma_vec,Particles,"MCMC",iters=10,ver=T)
        
        result<-c(K_eff/K,lZ_BH,lZ_Opt)
        return(result)
      },1:rep)
      #save means and sds
      listEstims<-c(listEstims,list(rbind(rowMeans(lZ_rep),rowSds(lZ_rep))))
      
      df_temp<-as.data.frame(t(lZ_rep)[,1])
      names(df_temp)<-"K_eff/K"
      p1<-ggplot(melt(df_temp), aes(x=variable, y=value)) + 
        geom_boxplot()
      p1<-p1 + geom_jitter(width = 0.05, height = 0.0) + ylim(0,1)+
              xlab("")+theme_bw()
      print(p1)
      dev.copy(pdf,paste0(plot_path,"K_eff_K",K,"_m",m,"_s",s,".pdf"),width=1, height=1.7)
      dev.off()
      
      df_temp<-as.data.frame(t(lZ_rep)[,-1])
      names(df_temp)<-c("logZ_BH","logZ_RB")
      p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) + 
        geom_boxplot()
      p2<- p2 + geom_jitter(width = 0.05, height = 0.0) + coord_cartesian(ylim=c(-1,1))+
        xlab("estimator") +theme_bw()
      print(p2)
      dev.copy(pdf,paste0(plot_path,"logZ_K",K,"_m",m,"_s",s,".pdf"),width=2, height=1.7)
      dev.off()
      
    }
  }
}

# Plots figure 3

K<-3e+4 #Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-3
m<-0.1; s<-2 #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0
gamma_vec<-(seq(0,1,l=21)^(5))[-1]

caseList<-list(c(0.5,2),c(0.2,20),c(0.2,Inf))

for(case in caseList){
  
  m<-case[1]
  s<-case[2]
  
  Particles<-r_q(N)
  plot(theta,exp(log_un_pi(theta)-lZ_true),type="l",col='blue',lwd=3,ylab='density')
  lines(density(Particles[,2]),col="red",lwd=3,lty=1)
  
  lZ_rep<-mcmapply(function(r){
    
    Particles<-r_q(N)
    Particles_alt<-altRep(Particles)
    
    lZ_BH<-BH_AIS_mod((seq(0,1,l=2)^(1))[-1],Particles,Particles_alt,"PG")
    lZ_mod3<-BH_AIS_mod(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=3)
    lZ_mod2<-BH_AIS_mod(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=2)
    lZ_mod1<-BH_AIS_mod(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=1)
    lZ_std_Gibbs<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG","Gibbs",iters=10)
    lZ_std_MCMC3<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=3)
    lZ_std_MCMC2<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=2)
    lZ_std_MCMC1<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=1)
    
    result<-c(lZ_BH,lZ_std_MCMC3,lZ_std_MCMC2,lZ_std_MCMC1,lZ_std_Gibbs,lZ_mod3,lZ_mod2,lZ_mod1)
    return(result)
  },1:rep)
  listEstims<-c(listEstims,list(rbind(rowMeans(lZ_rep),rowSds(lZ_rep))))
  
  df_temp<-as.data.frame(t(lZ_rep))
  names(df_temp)<-c("BH","AIS_M3","AIS_M2","AIS_M1","AIS_G","mAIS_M3","mAIS_M2","mAIS_M1")
  p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) + 
    geom_boxplot()
  p2<- p2 + geom_jitter(width = 0.05, height = 0.0) + 
    xlab("estimator")+theme_bw()+
    theme(axis.text.x = element_text(angle = 0))
  print(p2)
  dev.copy(pdf,paste0(plot_path,"BHAIS_logZ_K",K,"_m",m,"_s",s,".pdf"),width=6, height=1.7)
  dev.off()
}


# Plots comparison (sg) and (pg)

K<-3e+4 #Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-3
m<-0.1; s<-2 #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0
gamma_vec<-(seq(0,1,l=21)^(5))[-1]

caseList<-list(c(0.5,2),c(0.2,20),c(0.2,Inf))

for(case in caseList){
  
  m<-case[1]
  s<-case[2]
  
  Particles<-r_q(N)
  plot(theta,exp(log_un_pi(theta)-lZ_true),type="l",col='blue',lwd=3,ylab='density')
  lines(density(Particles[,2]),col="red",lwd=3,lty=1)
  
  lZ_rep<-mcmapply(function(r){
    
    Particles<-r_q(N)
    Particles_alt<-altRep(Particles)
    
    lZ_mod2_pg<-BH_AIS_mod(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=2)
    lZ_mod2_sg<-BH_AIS_mod(gamma_vec,Particles,Particles_alt,"SG","MCMC",iters=10,sigma_prop_MCMC=2)
    lZ_std_Gibbs_pg<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG","Gibbs",iters=10)
    lZ_std_Gibbs_sg<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"SG","Gibbs",iters=10)
    lZ_std_MCMC2_pg<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"PG","MCMC",iters=10,sigma_prop_MCMC=2)
    lZ_std_MCMC2_sg<-BH_AIS_std(gamma_vec,Particles,Particles_alt,"SG","MCMC",iters=10,sigma_prop_MCMC=2)
    
    result<-c(lZ_std_MCMC2_pg,lZ_std_MCMC2_sg,lZ_std_Gibbs_pg,lZ_std_Gibbs_sg,lZ_mod2_pg,lZ_mod2_sg)
    return(result)
  },1:rep)
  listEstims<-c(listEstims,list(rbind(rowMeans(lZ_rep),rowSds(lZ_rep))))
  
  df_temp<-as.data.frame(t(lZ_rep))
  names(df_temp)<-c("AIS_M_PG","AIS_M_SG","AIS_G_PG","AIS_G_SG","mAIS_M_PG","mAIS_M_SG")
  p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) + 
    geom_boxplot()
  p2<- p2 + geom_jitter(width = 0.05, height = 0.0) +theme_bw()+
    theme(axis.text.x = element_text(angle = 0)) +
    xlab("estimator")
  print(p2)
  dev.copy(pdf,paste0(plot_path,"BHAIS_PGSG_logZ_K",K,"_m",m,"_s",s,".pdf"),width=6.5, height=1.7)
  dev.off()
}


##### Intractable proposal

K<-5 #Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-3
m<-0.5; s<-2 #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0

Particles<-r_q(N)
plot(theta,exp(log_un_pi(theta)-lZ_true),type="l",col='blue',lwd=3,ylab='density',xlab="x")
lines(density(Particles[,2]),col="red",lwd=3,lty=1)
#dev.copy(pdf,paste0(plot_path,"density_K",K,"_m",m,"_s",s,".pdf"),width=4, height=5)
#dev.off()

caseList<-list(c(3),c(300))

for(case in caseList){

  K<- case[1]
  
  lZ_rep<-mcmapply(function(r){
    
    Particles<-r_q(N)
    lZ_ext<-IS_extended(Particles)
    
    Particles<-r_q(ceiling(N/K))
    lZ_opt<-IS_opt(Particles)
    
    result<-c(lZ_ext,lZ_opt)
    return(result)
  },1:rep)
  listEstims<-c(listEstims,list(rbind(rowMeans(lZ_rep),rowSds(lZ_rep))))
  
  df_temp<-as.data.frame(t(lZ_rep))
  names(df_temp)<-c("lZ_unif","lZ_opt")
  p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) +
    geom_boxplot()
  p2<-p2 + geom_jitter(width = 0.05, height = 0.0) +theme_bw()+
    theme(axis.text.x = element_text(angle = 0)) + coord_cartesian(ylim=c(-10,1)) +
    xlab("estimator")
  print(p2)
  dev.copy(pdf,paste0(plot_path,"IS_OptExt_comp_logZ_K",K,"_m",m,"_s",s,".pdf"),width=2, height=1.7)
  dev.off()
  
}


# Comparison of comb and opt

K<-300 #Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-3
m<-0.5; s<-2 #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0

N<-500

for(s in c(2,Inf))
for(N in c(500,3000))
for(K in c(30,300,500,3000,30000)){

  if(N==500){
    
    Particles<-r_q(N)
    
    df_temp<-data.frame(chain=rep("Proposal",N),
                        iteration=(1:N),
                        x=Particles[,2])
    
    p1<-ggplot(df_temp, aes(x,color=chain, fill=chain))+
      stat_function(data=data.frame(x = c(-25,25), chain=rep("Target",2)), mapping=aes(x), fun = function(x) exp(log_un_pi(x)-lZ_true), geom="area")+
      geom_area(stat="density",position="identity")+
      geom_line(stat="density", size=1)+
      scale_fill_manual(values=c(NA,"lightgray"))+
      scale_color_manual(values=c("black",NA))+
      coord_cartesian(xlim=c(-25,25), ylim=c(0,0.5))+
      ylab("density")+
      theme_classic()+
      theme(legend.title=element_blank(),legend.position="none")
    print(p1)
    dev.copy(pdf,paste0(plot_path,"density_K",K,"_m",m,"_s",s,".pdf"),width=2, height=1.7)
    dev.off()
  }
  
  lZ_rep<-mcmapply(function(r){
    
    Particles<-r_q(N)
    Particles_alt<-altRep(Particles)
    K_eff<-length(Particles_alt)
    
    lZ_BH<-BH_AIS_std(1,Particles,Particles_alt,"SG","None",iters=1)
    lZ_comb<-compute_lZ_G(Particles,Particles_alt)
    
    if(ceiling(N*K_eff/K)>1){
      
      Particles<-r_q(ceiling(N*K_eff/K))
      lZ_opt<-IS_opt(Particles)
    }
    else
      lZ_opt<--100
    
    result<-c(K_eff/K,K_eff/N,lZ_BH,lZ_opt,lZ_comb)
    return(result)
  },1:rep)
  listEstims<-c(listEstims,list(rbind(rowMeans(lZ_rep),rowSds(lZ_rep))))
  
  df_temp<-as.data.frame(t(lZ_rep)[,1:2])
  names(df_temp)<-c("K_eff/K","K_eff/N")
  p1<-ggplot(melt(df_temp), aes(x=variable, y=value)) +theme_bw()+
    geom_boxplot() + xlab("")
  p1<- p1 + geom_jitter(width = 0.05, height = 0.0) + ylim(0,1) +xlab("")
  print(p1)
  dev.copy(pdf,paste0(plot_path,"K_eff_N",N,"_K",K,"_m",m,"_s",s,".pdf"),width=1.7, height=1.7)
  dev.off()
  
  df_temp<-as.data.frame(t(lZ_rep)[,-(1:2)])
  names(df_temp)<-c("lZ_BH","lZ_opt","lZ_comb")
  p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) +theme_bw()+
    geom_boxplot()
  p2<- p2 + geom_jitter(width = 0.05, height = 0.0) + 
    theme(axis.text.x = element_text(angle = 0)) + coord_cartesian(ylim=c(min(df_temp,df_temp[which(df_temp[,2]==-100),2]),max(df_temp))) +
    xlab("estimator")
  print(p2)
  dev.copy(pdf,paste0(plot_path,"IS_Comb_comp_logZ_N",N,"_K",K,"_m",m,"_s",s,".pdf"),width=2.5, height=1.7)
  dev.off()
}

# Comparison of A1 and A2 

K<-3e+6#Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-3
m<-0.5; s<-2 #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0

N<-500

for(s in c(20,Inf))
for(N in c(500))
  for(K in c(3000,30000,3e+6)){
    
    if(s==20){
      
      Particles1<-r_q(N)
      s<-Inf
      Particles2<-r_q(N)
      s<-20
      
      df_temp1<-data.frame(chain=rep("Proposal1",N),
                          iteration=(1:N),
                          x=Particles1[,2])
      df_temp2<-data.frame(chain=rep("Proposal2",N),
                           iteration=(1:N),
                           x=Particles2[,2])
      df_temp<-rbind(df_temp1,df_temp2)
      
      p1<-ggplot(df_temp, aes(x,color=chain, fill=chain))+
        stat_function(data=data.frame(x = c(-30,30), chain=rep("Target",2)), mapping=aes(x), fun = function(x) exp(log_un_pi(x)-lZ_true), geom="area")+
        geom_area(stat="density",position="identity",linetype=3, size=0.9)+
        geom_line(aes(linetype=chain),stat="density", size=1)+
        scale_fill_manual(values=c(NA,NA,"lightgray"))+
        scale_color_manual(values=c("black","black",NA))+
        scale_linetype_manual(values = c("solid",NA,NA))+
        coord_cartesian(xlim=c(-25,25), ylim=c(0,0.5))+
        ylab("density")+
        theme_classic()+
        theme(legend.title=element_blank(),legend.position="none")
      print(p1)
      dev.copy(pdf,paste0(plot_path,"density2_K",K,"_m",m,"_s",s,".pdf"),width=2, height=1.7)
      dev.off()
      
      
    }
    
    lZ_rep<-mcmapply(function(r){
      
      Particles<-r_q(N)
      Particles_alt<-altRep(Particles)
      K_eff<-length(Particles_alt)
      
      lZ_BH<-BH_AIS_std(1,Particles,Particles_alt,"SG","None",iters=1)
      lZ_comb<-compute_lZ_G(Particles,Particles_alt)
      lZ_A1<-compute_lZ_A1(Particles,Particles_alt)[[1]]
      lZ_A2<-compute_lZ_A2(Particles,Particles_alt)[[1]]
      
      result<-c(K_eff/K,K_eff/N,lZ_BH,lZ_comb,lZ_A1,lZ_A2)
      return(result)
    },1:rep)
    df_temp<-as.data.frame(t(lZ_rep)[,-(1:2)])
    names(df_temp)<-c("lZ_BH","lZ_comb","lZ_GF1","lZ_GF2")
    p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) +theme_bw()+
      geom_boxplot()
    p2<- p2 + geom_jitter(width = 0.05, height = 0.0) + 
      theme(axis.text.x = element_text(angle = 0)) + coord_cartesian(ylim=c(min(df_temp[,-2],df_temp[which(df_temp[,2]!=-100),2]),max(df_temp)))+
      xlab("estimator")
    print(p2)
    dev.copy(pdf,paste0(plot_path,"IS_A1A2_comp_logZ_N",N,"_K",K,"_m",m,"_s",s,".pdf"),width=3, height=1.7)
    dev.off()
}


## AIS for GF

K<-3e+6 #Number of labels
mu_prop_min<--20; mu_prop_max<-20; sigma_prop<-sqrt(2); sigma_prop_MCMC<-3
m<-0.2; s<-Inf #Unifor: m=0.5, s=2, Binomial: m=0.5, s=Inf, Dicotomic: m=0.5, s=0
gamma_vec<-(seq(0,1,l=21)^(5))[-1]

N<-500
rep<-50

caseList<-list(c(0.3,50,3e+6),c(0.3,Inf,3e+6))
caseList<-list(c(0.3,20,3000),c(0.3,Inf,3000),
               c(0.3,20,30000),c(0.3,Inf,30000),
               c(0.3,20,3e+6),c(0.3,Inf,3e+6))

for(case in caseList){
  
  m<-case[1]
  s<-case[2]
  K<-case[3]

  if(s==Inf){
    
    m0<-m
    s0<-s
    
    Particles1<-r_q(N)
    s<-20
    Particles2<-r_q(N)
    m<-m0; s<-s0
    
    df_temp1<-data.frame(chain=rep("Proposal1",N),
                         iteration=(1:N),
                         x=Particles1[,2])
    df_temp2<-data.frame(chain=rep("Proposal2",N),
                         iteration=(1:N),
                         x=Particles2[,2])
    df_temp<-rbind(df_temp1,df_temp2)
    
    p1<-ggplot(df_temp, aes(x,color=chain, fill=chain))+
      stat_function(data=data.frame(x = c(-30,30), chain=rep("Target",2)), mapping=aes(x), fun = function(x) exp(log_un_pi(x)-lZ_true), geom="area")+
      geom_area(stat="density",position="identity",linetype=3, size=0.9)+
      geom_line(aes(linetype=chain),stat="density", size=1)+
      scale_fill_manual(values=c(NA,NA,"lightgray"))+
      scale_color_manual(values=c("black","black",NA))+
      scale_linetype_manual(values = c("solid",NA,NA))+
      coord_cartesian(xlim=c(-25,25), ylim=c(0,0.5))+
      ylab("density")+
      theme_classic()+
      theme(legend.title=element_blank(),legend.position="none")
    print(p1)
    dev.copy(pdf,paste0(plot_path,"density2_K",K,"_m",m,"_s",s,".pdf"),width=2, height=1.7)
    dev.off()
  }
  
  lZ_rep<-mcmapply(function(r){
    
    Particles<-r_q(N)
    Particles_alt<-altRep(Particles)
    
    lZ_BH<-BH_AIS_std(1,Particles,Particles_alt,"SG","None",iters=1)
    lZ_GF50<-GF_AIS_mod((seq(0,1,l=51+1)^(5))[-1],Particles,Particles_alt,"SG","MCMC",iters=10,sigma_prop_MCMC=3)
    lZ_GF20<-GF_AIS_mod((seq(0,1,l=21+1)^(5))[-1],Particles,Particles_alt,"SG","MCMC",iters=10,sigma_prop_MCMC=3)
    lZ_GF1<-GF_AIS_mod((seq(0,1,l=1+1)^(5))[-1],Particles,Particles_alt,"SG","MCMC",iters=10,sigma_prop_MCMC=3)
    
    result<-c(lZ_BH,lZ_GF1,lZ_GF20,lZ_GF50)
    return(result)
  },1:rep)
  
  df_temp<-as.data.frame(t(lZ_rep))
  names(df_temp)<-c("lZ_BH","GF_T1","GF_T21","GF_T51")
  p2<-ggplot(melt(df_temp), aes(x=variable, y=value)) + theme_bw()+
    geom_boxplot()
  p2<-p2 + geom_jitter(width = 0.05, height = 0.0) + theme(axis.text.x = element_text(angle = 0))+
    xlab("estimator")
  print(p2)
  dev.copy(pdf,paste0(plot_path,"GFAIS_SG_logZ_K",K,"_m",m,"_s",s,".pdf"),width=3, height=1.7)
  dev.off()

}

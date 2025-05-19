d<-5
alpha<-0.4
q<-qnorm(0.5*(1+(1-alpha)^(1/d)))
limSupp<-20

#Proposals parameters
C<-1e+6; k_pi<-1000
N_IS<-floor(C/(k_pi+1))

f<-function(K){
  
  N<-K
  return(N*k_pi+N*K-C)
}
plot(1:N_IS,f(1:N_IS),type="l")
rootInfo<-uniroot(f,lower=1,upper=N_IS)
K<-floor(rootInfo$root)
N_marg<-floor(C/(k_pi+K))

f<-function(N){
  
  return(N*k_pi+N*K*(1-(1-1/K)^N)-C)
}
plot(1:N_IS,f(1:N_IS),type="l")
rootInfo<-uniroot(f,lower=1,upper=N_IS)
N_BH<-floor(rootInfo$root)
c(N_marg,N_BH,N_IS)

C_marg<-N_marg*(k_pi+K)
C_BH<-N_BH*(k_pi+K*(1-(1-1/K)^N_BH))
C_IS<-N_IS*(k_pi+1)
c(C_marg,C_BH,C_IS)

print(detectCores())
options(mc.cores = detectCores()-1)
getOption("mc.cores",2L)

lUnDens_list<-list(ld1)
lWeights<-c(0)
lZ_true_vec<-c(ld1(0,"lZ",d=d))

lZ_true<-logSumExp(lZ_true_vec+lWeights)-logSumExp(lWeights)

#mu_list<-list(rep(0,d))
#Sigma_list<-list(diag(rep(10,d)))

nProps<-K
if(nProps==1){
  x1<-0
  x2<-.03*(100-x1^2)
}else {
  x1<-seq(-20,20,length.out=nProps)
  x2<-.03*(100-x1^2)
}

mu_list<-list()
Sigma_list<-list()
for(i in 1:nProps){
  
  mu_list[[i]]<-c(x1[i],x2[i],rep(0,d-2))
  Sigma_list[[i]]<-diag(c(10,1,rep(2,d-2)))
}
mu_IS<-list(rep(0,d)); Sigma_IS<-list(diag(c(100,10,rep(2,d-2))))

lZ_IS_Gauss<-logZ_BH2(N_IS,rep=50,mu_IS,Sigma_IS)
#lZ_IS_Unif<-logZ_IS(N_IS,rep=50,type="IS_unif")
#lZ_BH<-logZ_BH(N_BH,rep=30,K)
lZ_BH<-logZ_BH2(N_BH,rep=50,mu_list,Sigma_list)
lZ_marg<-logZ_marg2(N_marg,rep=50,mu_list,Sigma_list)
boxplot(cbind(lZ_IS_Gauss,lZ_BH,lZ_marg)); abline(h=lZ_true)

i1<-1; i2<-2
x1<-seq(-limSupp,limSupp,.1); x2<-x1
targetUnDens<-function(x,y){
  
  theta<-rep(0,d)
  theta[c(i1,i2)]<-c(x,y)
  
  return(logTargetUnDens(theta,lUnDens_list,lWeights))
}
z<-outer(x1,x1,FUN=Vectorize(targetUnDens))

plot(randProposal(N_BH,"marg",mu_list=mu_list,Sigma_list=Sigma_list)[[2]][,c(i1,i2)]); contour(x1,x2, z, drawlabels=FALSE, nlevels=100, add=T, col="blue")
plot(randProposal(N_IS,"marg",mu_list=mu_IS,Sigma_list=Sigma_IS)[[2]][,c(i1,i2)]); contour(x1,x2, z, drawlabels=FALSE, nlevels=100, add=T, col="blue")

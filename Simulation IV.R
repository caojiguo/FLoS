#=========== Simulation IV: the vairability of the estimators ========

#==== library pakages =====
# 
#install.packages("/path/SubsamplingFunPredictors_0.1.0.tar.gz", repos = NULL, type="source")
library(SubsamplingFunPredictors)
library(fda)
library(MASS)
library(mnormt) ## rmnorm: multinorm & rmt : multivariate t distribution
library(psych) ## the trace of matrix
library(wordspace) ## the rowNorms function
library(ggplot2)

# parameters setting
domain = c(0,1)
T = 101
nknots = 66
norder = 4
d = 3
n = 10^5
r = c(300,400,600,800,1000,1200,1400,1600,1800,2000,2200,2500,2800,3000)
r0 = r
K = ceiling(1.25*n^(0.25))
lambda = seq(0.1*n^(-3/8),1.5*n^(-3/8),length.out = 10)

# scenario
aind = 1
a = 0
b = 6
df = NA
family = "Binomial"

# generate X
knots    = seq(domain[1],domain[2], length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)

tobs = seq(domain[1],domain[2],length.out = T)
basismat = eval.basis(tobs, basis)

X = a_fun(n,nbasis,aind,a,b,df) %*% t(basismat)
matplot(tobs, t(X[sample(1:n,3,replace = FALSE),]),type="l",xlab="t",ylab="x(t)")

# new design matrix and smoothness matrix
NV = compute.NV(X,K,d,domain)
N = NV$N
V = NV$V
N_norm = NV$N_norm
basismat = NV$basismat


# functional coefficient
betaeval = 8*sin(0.85*pi*tobs)


# y0 the signals
h   = (domain[2]-domain[1])/(T-1)
cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
y0  = h/3*X%*%diag(cef)%*%betaeval
boxplot(y0)


# generate Y
prob = psi(y0,family)
Y = generator_Y_FGLM(n,prob,family)
table(Y)/10^5


# IMSE result
S = 300
IMSE_1 = array(0,dim = c(S,2,length(r)))
PCC_1 = array(0,dim = c(S,2,length(r)))


for (i in 1:S){
   for (j in 1:length(r)){
    try({
      # BIC
      FLoS_result = FLoS_FGLM_BIC(N,N_norm,Y,r[j],r0[j],lambda,V,family)
      lambda_FLoS = FLoS_result$lambda_FLoS
      c_FLoS = FLoS_result$c_FLoS
      c0 = FLoS_result$c0
      c_uni = Unisub_FGLM(N, Y, r[j], lambda_FLoS, V,c0,family)
      
      
      # estimate beta(t)
      beta_true = 8*sin(0.85*pi*tobs)
      beta_FLoS = basismat%*%c_FLoS
      
      beta_uni = basismat%*%c_uni
      
      
      
      IMSE_FLoS = sqrt(mean((beta_FLoS-beta_true)^2))
      
      IMSE_uni = sqrt(mean((beta_uni-beta_true)^2))
      
      
      # Proportions of correct classification (PCC)
      if(family == "Binomial"){
        Y_prob_FLoS <- 1/(1 + exp(-N %*% c_FLoS))
        Y_hat_FLoS <- 1*(Y_prob_FLoS>0.5)
        PCC_FLoS = sum(Y==Y_hat_FLoS )/n
        Y_prob_unif <- 1/(1 + exp(-N %*% c_uni))
        Y_hat_unif <- 1*(Y_prob_unif>0.5)
        PCC_unif = sum(Y==Y_hat_unif )/n
        
        result = list(IMSE = c(IMSE_FLoS,  IMSE_uni),PCC = c(PCC_FLoS,PCC_unif))
      } else { result = list(IMSE = c(IMSE_FLoS,  IMSE_uni))}
      IMSE_1[i,,j] = result$IMSE
      PCC_1[i,,j] = result$PCC
      print(c(i,j))
    }, silent=TRUE)
   }
}

#===== plots =====
# IMSE

Method = c(rep("FLoS",length(r)),rep("UNIS",length(r)))
sIMSE_1 = apply(IMSE_1,c(2,3),sd)
datag = data.frame(r = rep(r,2), Method = Method, value = as.vector(t(sIMSE_1)))
ggplot(data=datag[-c(1,15),], aes(x=r, y=value, group=Method,colour=Method)) +
  geom_point(aes(shape = Method))+
  geom_line(aes(linetype = Method))+
  scale_color_manual(values = c("red", "black"))+
  labs(x="subsample size", y="Sd(RIMSE)")+
  theme(legend.position=c(0.8,0.8))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



# boxplot
datanew = data.frame(FLoSIMSE =as.vector(IMSE_1[,1,]),r=rep(as.factor(r),each = 300))
boxplot(IMSE_1[,1,14])

boxplot(FLoSIMSE~r,data= datanew,ylab = "RIMSE",xlab = "subsample size ")
ggplot(datanew, aes(x = r, y = FLoSIMSE, fill = r)) + 
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")+
  labs(x="subsample size", y="RIMSE")+
  theme(legend.position=c(1,5))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(datanew, aes(x = r, y = FLoSIMSE, fill = r)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2)+
  labs(x="subsample size", y="RIMSE")+
  theme(legend.position=c(1,5))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



# beta result
r = c(1000,2000,2500,3000)
r0 = r
S = 100
beta_result = array(0,dim = c(S,T,length(r)))

for (i in 1:S){
  for (j in 1:length(r)){
    try({
      # BIC
      FLoS_result = FLoS_FGLM_BIC(N,N_norm,Y,r[j],r0[j],lambda,V,family)
      lambda_FLoS = FLoS_result$lambda_FLoS
      c_FLoS = FLoS_result$c_FLoS
      c0 = FLoS_result$c0
      c_uni = Unisub_FGLM(N, Y, r[j], lambda_FLoS, V,c0,family)
      
      
      # estimate beta(t)
      beta_true = 8*sin(0.85*pi*tobs)
      beta_FLoS = basismat%*%c_FLoS
      
      beta_result[i,,j] = beta_FLoS
      print(c(i,j))
    }, silent=TRUE)
  }
}


# plot of 100 lines
beta_1000 = beta_result[,,1]
beta_2000 = beta_result[,,2]
beta_2500 = beta_result[,,3]
beta_3000 = beta_result[,,4]

beta_full = data.frame(tobs = tobs,beta=betaeval,trial = rep(paste("trial",0),101))
data_beta_1000 = data.frame(tobs = rep(tobs,100),beta = as.vector(t(beta_1000)),
                            trial = rep(paste("trial",1:100),each=101))
data_beta_2000 = data.frame(tobs = rep(tobs,100),beta = as.vector(t(beta_2000)),
                            trial = rep(paste("trial",1:100),each=101))
data_beta_2500 = data.frame(tobs = rep(tobs,100),beta = as.vector(t(beta_2500)),
                            trial = rep(paste("trial",1:100),each=101))

data_beta_3000 = data.frame(tobs = rep(tobs,100),beta = as.vector(t(beta_3000)),
                            trial = rep(paste("trial",1:100),each=101))



ggplot(data=data_beta_1000, aes(x=tobs, y=beta, group=trial,colour=trial)) +
  geom_line(linetype =3)+
  scale_colour_grey()+
  geom_line(data=beta_full, aes(x=tobs, y=beta),size = 1,linetype = 1,color = "red")+
  scale_y_continuous(limits = c(-1.5,9))+
  labs(x="t", y=expression(beta(t)))+
  theme(legend.position=c(2,2))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(data=data_beta_2000, aes(x=tobs, y=beta, group=trial,colour=trial)) +
  geom_line(linetype =3)+
  scale_colour_grey()+
  geom_line(data=beta_full, aes(x=tobs, y=beta),size = 1,linetype = 1,color = "red")+
  scale_y_continuous(limits = c(-1.5,9))+
  labs(x="t", y=expression(beta(t)))+
  theme(legend.position=c(2,2))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(data=data_beta_2500, aes(x=tobs, y=beta, group=trial,colour=trial)) +
  geom_line(linetype =3)+
  scale_colour_grey()+
  geom_line(data=beta_full, aes(x=tobs, y=beta),size = 1,linetype = 1,color = "red")+
  scale_y_continuous(limits = c(-1.5, 9))+
  labs(x="t", y=expression(beta(t)))+
  theme(legend.position=c(2,2))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(data=data_beta_3000, aes(x=tobs, y=beta, group=trial,colour=trial)) +
  geom_line(linetype =3)+
  scale_colour_grey()+
  geom_line(data=beta_full, aes(x=tobs, y=beta),size = 1,linetype = 1,color = "red")+
  scale_y_continuous(limits = c(-1.5, 9))+
  labs(x="t", y=expression(beta(t)))+
  theme(legend.position=c(2,2))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))




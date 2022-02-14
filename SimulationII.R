#==== Simulation III: Scenario I & n= 1000000 ========

#==== library pakages =====
# install.packages("SubsamplingFunPredictors_0.1.0.tar.gz", repos = NULL, type="source")
library(SubsamplingFunPredictors)
library(fda)
library(MASS)
library(mnormt) ## rmnorm: multinorm & rmt : multivariate t distribution
library(psych) ## the trace of matrix
library(wordspace) ## the rowNorms function

# parameters setting
domain = c(0,1)
T = 101
d = 3
n = 10^6
r = 1000*(1:6)
r0 = r/2
K = ceiling(1.25*n^(0.25))
lambda = seq(0.1*n^(-3/8),1.5*n^(-3/8),length.out = 10)

# scenario
aind = 1
a = 0
b = 1
df = NA
family = "Possion"
# family = "Binomial"


# IMSE result
IMSE_1 = array(0,dim = c(100,2,6))
PCC_1 = array(0,dim = c(100,2,6))
for (i in 1:100){
  data = FGLMR.data.generator.bsplines(n=n,nknots=66,norder=4,T = T, domain=c(0,1), aind=aind,a=a,b=b,df=df)
  X   = data$X
  y0   = data$y0

  tobs   = seq(domain[1],domain[2],length.out = T)

  # generate Y
  prob = psi(y0,family)
  Y = generator_Y_FGLM(n,prob,family)


  # new design matrix and smoothness matrix
  NV = compute.NV(X,K,d,domain)
  N = NV$N
  V = NV$V
  N_norm = NV$N_norm
  basismat = NV$basismat

  for (j in 1:6){
    try({
      # BIC
      FLoS_result = FLoS_FGLM_BIC(N,N_norm,Y,r[j],r0[j],lambda,V,family)
      lambda_FLoS = FLoS_result$lambda_FLoS
      c_FLoS = FLoS_result$c_FLoS
      c0 = FLoS_result$c0
      c_uni = Unisub_FGLM(N, Y, r[j], lambda_FLoS, V,c0,family)


      # estimate beta(t)
      beta_true = sin(0.5*pi*tobs)

      beta_FLoS = basismat%*%c_FLoS

      beta_uni = basismat%*%c_uni


      IMSE_FLoS = sqrt(mean((beta_FLoS-beta_true)^2))

      IMSE_uni = sqrt(mean((beta_uni-beta_true)^2))

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
      # PCC_1[i,,j] = result$PCC
      print(c(i,j))
    }, silent=TRUE)
  }
}

#===== plots =====
# IMSE
Method = c(rep("FLoS",6),rep("UNIS",6))
mIMSE_1 = apply(IMSE_1,c(2,3),mean)
datag = data.frame(r = rep(r,2), Method = Method, value = as.vector(t(mIMSE_1)))
ggplot(data=datag, aes(x=r, y=value, group=Method,colour=Method)) +
  geom_point(aes(shape = Method))+
  geom_line(aes(linetype = Method))+
  scale_color_manual(values = c("black", "red"))+
  # ggtitle("Scenario I")+
  labs(x="subsample size", y="IMSE")+
  theme(legend.position=c(0.8,0.8))+
  theme(legend.key = element_blank())+
  theme(legend.background = element_blank())+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))





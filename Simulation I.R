#=========== Simulation I: Scenario I & n= 100000 ========

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
r = c(300,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
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
matplot(tobs, t(X[sample(1:n,10,replace = FALSE),]),type="l",xlab="t",ylab="x(t)")

# new design matrix and smoothness matrix
NV = compute.NV(X,K,d,domain)
N = NV$N
V = NV$V
N_norm = NV$N_norm
basismat = NV$basismat


# functional coefficient
betaeval = 8*sin(0.85*pi*tobs)
plot(tobs,betaeval)


# y0 the signals
h   = (domain[2]-domain[1])/(T-1)
cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
y0  = h/3*X%*%diag(cef)%*%betaeval


# IMSE result
S = 300
IMSE_1 = array(0,dim = c(S,2,length(r)))
PCC_1 = array(0,dim = c(S,2,length(r)))

for (i in 1:S){
  # generate Y
  prob = psi(y0,family)
  Y = generator_Y_FGLM(n,prob,family)
  
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

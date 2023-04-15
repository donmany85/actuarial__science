#priory rating for cross sectioal data

# Experience Rating for panel data

#from panel to multilevel data



#Linear odels for longitudinal data



library(CASdatasets)

data(usmassBI) 
data(usmassBI2)

AutoClaim=usmassBI2

names(AutoClaim)



AutoClaim[1:12,]



AutoClaimIn=subset(AutoClaim, YEAR<1998)

#Multiple time series plot



plot(AC~YEAR, data=AutoClaimIn, ylab='Average Claim', xlab='Year')
for (i in AutoClaimIn$TOWNCODE){
  lines(AC~YEAR, data=subset(AutoClaimIn, TOWNCODE==i))}



#scatter plot to Explore relations
AutoClaimIn$lnPCI=log(AutoClaimIn$PCI)
AutoClaimIn$lnPPSM=log(AutoClaimIn$PPSM)

plot(AC~lnPCI, data=AutoClaimIn, ylab='Average Claim', xlab='PCI')
for (i in AutoClaimIn$TOWNCODE){
  lines(AC~lnPCI, data=subset(AutoClaimIn, TOWNCODE==i))}


plot(AC~lnPPSM, data=AutoClaimIn, ylab='Average Claim', xlab='PPSM')
for (i in AutoClaimIn$TOWNCODE){
  lines(AC~lnPPSM, data=subset(AutoClaimIn, TOWNCODE==i))}


AutoClaimIn$YEAR=AutoClaimIn$YEAR-1992
Pool.fit=lm(AC~lnPCI+lnPCI+lnPPSM+YEAR, data=AutoClaimIn)

summary(Pool.fit)



#Fixed Effects Models

#Basic fixed-effects model

FE.fit=lm(AC~factor(TOWNCODE)+lnPCI+lnPPSM+YEAR-1, data=AutoClaimIn)
summary(FE.fit)


anova(Pool.fit,FE.fit)


# Models with Serial Correlation

AutoClaimIn$rPool=resid(Pool.fit)
rvec=cbind(subset(AutoClaimIn,YEAR==1)$rPool, subset(AutoClaimIn, YEAR==2)$rPool,
           subset(AutoClaimIn, YEAR==3)$rPool, subset(AutoClaimIn, YEAR==4)$rPool,
           subset(AutoClaimIn, YEAR==5)$rPool)
cor(rvec)

library(nlme)
SCex.fit=gls(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
             correlation=corCompSymm(form=~1|TOWNCODE))
summary(SCex.fit)

intervals(SCex.fit, which='var-cov')
getVarCov(SCex.fit)


#AR(1)
SCar.fit=gls(AC~lnPCI+lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
             correlation=corAR1(form=~1|TOWNCODE))
summary(SCar.fit)

intervals(SCar.fit, which='var-cov')
getVarCov(SCar.fit)


#Unstructured
SCun.fit=gls(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
             correlation  =corSymm(form=~1|TOWNCODE))
summary(SCun.fit)
intervals(SCun.fit, which='var-cov')
getVarCov(SCun.fit)

#Likelihood ratio test
SCex.fit.ml=gls(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
                correlation=corCompSymm(form=~1|TOWNCODE), method='ML')
SCar.fit.ml=gls(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
                correlation = corAR1(form=~1|TOWNCODE), method='ML')
SCun.fit.ml=gls(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
                correlation=corSymm(form=~1|TOWNCODE), method='ML')

anova(SCex.fit.ml,Pool.fit)

anova(SCar.fit.ml, Pool.fit)
anova(SCun.fit.ml, Pool.fit)



# library(nlme)

#error -components model

EC.fit=lme(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn, random=~1|TOWNCODE)

summary(EC.fit)

intervals(EC.fit, which='var-cov')


#Pooling test

tcode=unique(AutoClaimIn$TOWNCODE)
n=length(tcode)
N=nrow(AutoClaimIn)
T=rep(NA,n)
s=rep(NA,n)
for (i in 1:n){
  T[i]=nrow(subset(AutoClaimIn, TOWNCODE==tcode[i]))
  s[i]=(sum(subset(AutoClaimIn, TOWNCODE==tcode[i])$rPool)^2
        -sum(subset(AutoClaimIn, TOWNCODE==tcode[i])$rPool^2))/T[i]/(T[i]-1)}
TS=(sum(s*sqrt(T*(T-1)))*N/sum(AutoClaimIn$rPool^2))^2/2/n
TS


##Error component with AR1
RE.fit=update(EC.fit, correlation=corAR1(form=~1|TOWNCODE))
summary(RE.fit)


intervals(RE.fit, which='var-cov')


##Get variance components

getVarCov(RE.fit)

getVarCov(RE.fit, type='conditional')


getVarCov(RE.fit, type='marginal')


#likelihood ratio test

EC.fit.ml=lme(AC~lnPCI+lnPPSM+YEAR, data=AutoClaimIn,
              random=~1|TOWNCODE, method='ML')
RE.fit.ml=update(EC.fit, correlation=corAR1(form=~1|TOWNCODE), method='ML')
anova(EC.fit.ml, RE.fit.ml)

##Hausman test

Var.FE=vcov(FE.fit)[-(1:n),-(1:n)]
Var.EC=vcov(EC.fit)[-1,-1]
beta.FE=coef(FE.fit)[-(1:n)]
beta.EC=fixef(EC.fit)[-1]
ChiSq=t(beta.FE-beta.EC)%*%solve(Var.FE-Var.EC)%*%(beta.FE-beta.EC)
ChiSq



#prediction


#BLUP

alpha.BLUP=ranef(EC.fit)
beta.GLS=fixef(EC.fit)
resid.BLUP=residuals(EC.fit, type='response')
rstandard.BLUP=residuals(EC.fit, type='normalized')
alpha.BLUP


#Use data of year 1998 for validation
AutoClaimOut=subset(AutoClaim, YEAR==1998)

#Define new variables
AutoClaimOut$lnPCI=log(AutoClaimOut$PCI)
AutoClaimOut$lnPPSM=log(AutoClaimOut$PPSM)
AutoClaimOut$YEAR=AutoClaimOut$YEAR-1992

#Compare models Pool.fit, SCar.fit, FE.fit, Ec.fit, Re.fit and Fear.fit
#Fixed -effects model with AR(1)

FEar.fit =gls(AC~factor(TOWNCODE)+lnPCI+lnPPSM+YEAR-1,
              data=AutoClaimIn, correlation=corAR1(form=~1|TOWNCODE))
FEar.fit.ml=gls(AC~factor(TOWNCODE)+lnPCI+lnPPSM+YEAR-1,
                data=AutoClaimIn,correlation=corAR1(form=~1|TOWNCODE), method='ML')


#prediction

Xmat=cbind(rep(1,nrow(AutoClaimOut)), AutoClaimOut$lnPCI,
           AutoClaimOut$lnPPSM, AutoClaimOut$YEAR)
beta.Pool=coef(Pool.fit)
pred.Pool=Xmat%*%beta.Pool
MSPE.Pool=sum((pred.Pool-AutoClaimOut$AC)^2)
MAPE.Pool=sum(abs(pred.Pool-AutoClaimOut$AC))


beta.SCar=coef(SCar.fit)
pred.SCar=Xmat%*%beta.SCar
MSPE.SCar=sum((pred.SCar-AutoClaimOut$AC)^2)
MAPE.SCar=sum(abs(pred.SCar-AutoClaimOut$AC))

beta.FE=coef(FEar.fit)[-(1:29)]
pred.FE=coef(FE.fit)[1:29]+Xmat[,-1]%*%beta.FE
MSPE.FE=sum((pred.FE-AutoClaimOut$AC)^2)
MAPE.FE=sum(abs(pred.FE-AutoClaimOut$AC))

beta.FEar=coef(FEar.fit)[-(1:29)]
pred.FEar=coef(FEar.fit)[1:29]+Xmat[,-1]%*%beta.FEar
MSPE.FEar=sum((pred.FEar-AutoClaimOut$AC)^2)
MAPE.FEar=sum(abs(pred.FEar-AutoClaimOut$AC))

alpha.EC=ranef(EC.fit)

beta.EC=fixef(EC.fit)
pred.EC=alpha.EC+Xmat%*%beta.EC
MSPE.EC=sum((pred.EC-AutoClaimOut$AC)^2)
MAPE.EC=sum(abs(pred.EC-AutoClaimOut$AC))

alpha.RE=ranef(RE.fit)
beta.RE=fixef(RE.fit)
pred.RE=alpha.RE+Xmat%*%beta.RE
MSPE.RE=sum((pred.RE-AutoClaimOut$AC)^2)
MAPE.RE=sum(abs(pred.RE-AutoClaimOut$AC))




##Generalized Linear Models for Longitudinal Data

#Specifying Generalized Linear Models with Random Effects



#experience rating with bonus-malus scales in R

Pmatrix=function(th){
  P=matrix(nrow=6,ncol=6,data=0)
  P[1,1]=P[2,1]=P[3,2]=P[4,3]=P[5,4]=P[6,5]=exp(-th)
  P[,6]=1-exp(-th)
  return(P)
}

lim.distr=
  function(matrix){
    et=matrix(nrow=1, ncol=dim(matrix)[2], data=1)
    E=matrix(nrow=dim(matrix)[1], ncol=dim(matrix)[2],data=1)
    mat=diag(dim(matrix)[1])-matrix+E
    inverse.mat=solve(mat)
    p=et%*%inverse.mat
    return(p)  }

P=Pmatrix(0.1)
P
pi=lim.distr(P)
pi



## without a priori ratemaking

a.hat= 0.8888

lambda.hat=0.1474
int1=
  function(theta, s, a, lambda){
    a=a.hat
    lambda=lambda.hat
    f.dist=gamma(a)^(-1)*a^a*theta^(a-1)*exp(-a*theta)
    p=lim.distr(Pmatrix((lambda*theta)))
    return(theta*p[1,s+1]*f.dist)  }
P1=matrix(nrow=1, ncol=6, data=0)
for(i in 0:5) P1[1,i+1]=integrate(Vectorize(int1), lower=0, upper=Inf, s=i)$value

int2=
  function(theta,s,a,lambda){
    a=a.hat
    lambda=lambda.hat
    f.dist=gamma(a)^(-1)*a^a*theta^(a-1)*exp(-a*theta)
    P=lim.distr(Pmatrix((lambda*theta)))
    return(P[1,s+1]*f.dist)  }

P2=matrix(nrow=1, ncol=6,data=0)
for (i in 0:5) P2[1,i+1]=integrate(Vectorize(int2), lower=0, upper=Inf, s=i)$value

R=P1/P2
R



##nagative binomial regression

lambda= c(0.1176,0.1408,0.1897,0.2272,0.1457,0.1746,0.2351,0.2816,
          0.1761,0.2109,0.2840,0.3402,0.2182,0.2614,0.3520,0.0928,
          0.1112,0.1498,0.1794,0.1151,0.1378,0.1856,0.2223)
weights = c(0.1049,0.1396,0.0398,0.0705,0.0076,0.0122,0.0013,0.0014,
            0.0293,0.0299,0.0152,0.0242,0.0007,0.0009,0.0002,0.1338,
            0.1973,0.0294,0.0661,0.0372,0.0517,0.0025,0.0044)

a=1.065
n=length(weights)

int3=
  function(theta,lambda, a, l){
    P=lim.distr(Pmatrix((lambda*theta)))
    f.dist=gamma(a)^(-1)*a^a*theta^(a-1)*exp(-a*theta)
    return(theta*P[1,l+1]*f.dist)  }

int4=
  function(theta,lambda, a, l){
    P=lim.distr(Pmatrix((lambda*theta)))
    f.dist=gamma(a)^(-1)*a^a*theta^(a-1)*exp(-a*theta)
    return(P[1,l+1]*f.dist)  }

teller1=teller2=noemer=array(dim=6,data=0)
result1=result2=array(dim=6,data=0)

for (i in 0:5){
  b=c=array(dim=n, data=0)
  for (j in 1:n){
    b[j]=integrate(Vectorize(int3), lower=0, upper=Inf, lambda=lambda[j], a=a,l=i)$value
    c[j]=integrate(Vectorize(int4), lower=0, upper=Inf, lambda=lambda[j], a=a,l=i)$value  }
  teller1[i+1]=b%*% weights
  noemer[i+1]=c%*%weights
  R=teller1/noemer}

R
#Probability Distributions in Actuarial Science

# Pearson system
install.packages('PearsonDS')
library(PearsonDS)

x=seq(-1,6, by=1e-3)
y0=dpearson0(x,2,1/2)
y1=dpearsonI(x,1.5,2,0,2)
y2=dpearsonII(x,2,0,1)
y3=dpearsonIII(x,3,0,1/2)
y4=dpearsonIV(x,2.5,1/3,1,2/3)
y5=dpearsonV(x,2.5,-1,1)
y6=dpearsonVI(x,1/2,2/3,2,1)
y7=dpearsonVII(x,3,4,1/2)
plot(x,y0,type='l', ylim=range(y0, y1, y2, y3,y4,y5,y7), 
     ylab='f(x)', main='The Pearson distribution system')
lines(x[y1 !=0], y1[y1!=0], lty=2)
lines(x[y2 !=0], y2[y2!=0], lty=3)
lines(x[y3 !=0], y3[y3!=0], lty=4)
lines(x,y4, col='grey')
lines(x,y5, col='grey', lty=2)
lines(x[y6 !=0], y6[y6!=0], lty=3)
lines(x[y7 !=0], y7[y7!=0], lty=4)
legend('topright', leg=paste('Pearson', 0:7), lty=1:4,
       col=c(rep('black',4), rep('grey',4)))



#Exponential family


dgamma(1:2, shape=2, rate=3/2)

pgamma(1:2, shape=2,rate=3/2)
qgamma(1/2, shape=2, rate=3/2)

set.seed(1)
rgamma(5,shape=2, rate=3/2)


## cran.r-project.org/web/views/Distributions.html
### acuar, actuDistns
##actuar provides the raw moment E(X^k) the limited expected values E(min(X,l)k)
# and the moment generating functions E(e^tX) for manydistributions in three
#dedicated functions mfoo, levfoo mgffoo

##when  on a particular problem all classical distributions have been exhausted
#it is some times appropriate to create new probability distributions.

#typical transformations of a random variable X are listed:

##(i) Translation X+c(e.g the shifted lognormal distribution)
##(ii) Scaling lambda*X
##(iii) PoWer X^alpha(e.g the generalized beta type 1distribution)
##(iv) Inverse 1/X (e.g. the inverse gamma distribution)
##(v) The logarithm log(X) (e.g. the loglogistic distribution)
##(vi) Exponential exp(X) and
##(vii) the odd ratio X/(1-X)  (e.g. the beta type 2 distribution)


f = function(x) dgamma(x,2)
f1= function(x) f(x-1)
f2= function(x) f(x/2)/2
f3= function(x) 2*x*f(x^2)
f4= function(x) f(1/x)/x^2
f5= function(x) f(exp(x))*exp(x)
f6= function(x) f(log(x))/x

x=seq(0,10,by=.025)
plot(x,f(x), ylim=c(0,1.3), xlim=c(0,10), main='Theoretial densities',
     lwd=2, type='l', xlab='x', ylab='')
lines(x,f1(x),lty=2, lwd=2)
lines(x,f2(x), lty=3, lwd=2)
lines(x, f3(x), lty=4, lwd=2)
lines(x,f4(x), lty=1, col='grey', lwd=2)
lines(x,f5(x), lty=2, col='grey', lwd=2)
lines(x,f6(x), lty=3, col='grey', lwd=2)
legend('topright', lty=1:4, col=c(rep('black',4), rep('grey', 3)),
       leg=c('X','X+1','2X','sqrt(X)','1/X', 'log(X)', 'exp(X)'))

##simulations and visualize kernel-based densities

set.seed(123)
x= rgamma(100,2)
x1= x+1
x2= 2*x
x3= sqrt(x)
x4= 1/x
x5= log(x)
x6= exp(x)
plot(density(x), ylim=c(0,1), xlim=c(0,10), main='Empirical densities',
     lwd=2, xlab='x', ylab='f_X(x)')
lines(density(x1), lty=2, lwd=2)
lines(density(x2), lty=3, lwd=2)
lines(density(x3), lty=4, lwd=2)
lines(density(x4), lty=1, col='grey', lwd=2)
lines(density(x5), lty=2, col='grey', lwd=2)
lines(density(x6), lty=3, col='grey', lwd=2)

dpois(0:2, lambda=3)

ppois(1:2, lambda=3)

qpois(1/2, lambda=3)

rpois(5,lambda=3)


##zero-modified poisson

dpoisZM=function(x,prob,lambda){
  prob*(x==0)+(1-prob)*(x>0)*dpois(x-1,lambda)}
ppoisZM=function(q,prob,lambda){
  prob*(q>=0)+(1-prob)*(q>0)*ppois(q-1,lambda)}
qpoisZM=function(p,prob,lambda){
  ifelse(p<=prob,0,1+qpois((p-prob)/(1-prob),lambda))}
rpoisZM=function(n,prob,lambda){
  (1-rbinom(n,1,prob))*(rpois(n,lambda)+1)}




## Mixed type distribution

# one-modified betta
dbeta0M=function(x,prob,a,b)
{dbeta(x,a,b)*(1-prob)*(x!=1)+prob*(x==1)}
pbeta0M=function(q,prob,a,b)
{pbeta(q,a,b)*(1-prob)+prob*(q>=1)}

##MBBEFD distribution(some context of reinsurance treaties)
dMBBEFD=function(x,a,b)
{-a*(a+1)*b^x*log(b)/(a+b^x)^2+(a+1)*b/(a+b)*(x==1)}
pMBBEFD=function(x,a,b)
{a*((a+1)/(a+b^x)-1)*(x<1)+1*(x>=1)}

# install.packages('actuar')
library(actuar)
dmixgampar=function(x,prob,nu,lambda,alpha,theta){
  prob*dgamma(x,nu,lambda)+(1-prob)*dpareto(x,alpha,theta)}
pmixgampar=function(q,prob,nu,lambda,alpha,theta){
  prob*pgamma(q,nu,lambda)+(1-prob)*ppareto(q,alpha,theta)}

M=matrix(0,3,3)
diag(M)=-2
diag(M[1:(nrow(M)-1),2:ncol(M)])=-2
M

set.seed(123)
rphtype(5,prob=c(1,0,0), rates=M)

# install.packages('distr')
# install.packages('distrEx')
library(distr)
library(distrEx)
X=Norm(mean=5,sd=2)
X

#rstudio에서 q를 쓰면 quit이라 q.l사용
q.l(X)
function(p, lower.tail=TRUE, log.p=FALSE)
{qnorm(p,mean=5, sd=2, lower.tail=lower.tail, log.p=log.p)}
q.l(X)(0.25)
mean(X)



N=DiscreteDistribution(supp=c(1,2,4,9), prob=c(.2,.4,.35,.05))

plot(N)

X1=Norm(mean=5,sd=2)
X2=Norm(mean=2,sd=1)


S=X1+X2
plot(S)
U=Unif(Min=0, Max=1)
N=DiscreteDistribution(supp=c(1,2,4,9), prob=c(.2,.4,.35,.05))
Z=U^N
plot(Z)



N=DiscreteDistribution(supp=c(0,1,2,4), prob=c(.2,.4,.35,.05))
Z=U^N
Z
plot(Z)

CP=CompoundDistribution(Pois(),Gammad())
CP

data(itamtplcost)
# install.packages('fitdistrplus')
library(MASS)
library(survival)
library(fitdistrplus)

x=itamtplcost$UltimateCost/10^6
summary(x)
fgamMLE=fitdist(x,'gamma', method='mle')
fgamMLE


summary(fgamMLE)



##Moment Matching Estimation
fgamMME=fitdist(x,'gamma', method='mme')
cbind(MLE=fgamMLE$estimate, MME=fgamMME$estimate)


##Quantile Matching Estimation
fgamQME=fitdist(x,'gamma', method='qme', probs=c(1/3,2/3))
cbind(MLE=fgamMLE$estimate, MME=fgamMME$estimate, QME=fgamQME$estimate)


#Maximum Goodness-of Fit Estimation
fgamMGE=fitdist(x,'gamma', method='mge', gof='CvM')
cbind(MLE=fgamMLE$estimate, MME=fgamMME$estimate,
      QME=fgamQME$estimate, MGE=fgamMGE$estimate)


#Measures of Adequacy

##Histogram and Empirical Densities

txt=c('MLE', 'MME', 'QME(1/3,2/3)', 'MGE-CvM')
par(mfrow=c(1,1))
denscomp(list(fgamMLE, fgamMME, fgamQME, fgamMGE), legendtext=txt,
         fitcol='black', main='Histogram and fitted gamma densities')

hist(x,prob=TRUE, ylim=c(0,1))
lines(density(x), lty=5)


#Distribution Function Plot

data(danishuni)
x=danishuni$Loss
fgam=fitdist(x,'gamma', lower=0)
fpar=fitdist(x,'pareto', start=list(shape=2, scale=2), lower=0)
fmixgampar=fitdist(x,'mixgampar', start=
                     list(prob=1/2, nu=1, lambda=1, alpha=2, theta=2), lower=0)
cbind(SINGLE=c(NA, fgam$estimate, fpar$estimate),
      MIXTURE=fmixgampar$estimate)

fburr=fitdist(x,'burr', start=list(shape1=2, shape2=2,
                                   scale=2), lower=c(0.1,1/2,0))
fburr$estimate

cdfcomp(list(fgam, fpar, fmixgampar, fburr), xlogscale=TRUE,
        datapch='.', datacol='grey', fitcol='black', fitlty=2:5,
        legendtext=c('gamma', 'Pareto', 'Par-gam', 'Burr'),
        main='Fitted CDFs on danish')


##QQplot PPplot
qmixgampar=function(p,prob,nu,lambda, alpha, theta)
{
  L2=function(q,p)
    (p-pmixgampar(q,prob, nu, lambda, alpha, theta))^2
  sapply(p, function(p) optimize(L2, c(0,10^3), p=p)$minimum)
}


ppcomp(list(fgam, fpar, fmixgampar, fburr), xlogscale=TRUE,
       ylogscale=TRUE, fitcol='black', main='PP-plot on danish',
       legendtext=c('gamma', 'Pareto', 'Par-gam', 'Burr'), fitpch=1:4)

qqcomp(list(fgam, fpar, fmixgampar, fburr), xlogscale=TRUE,
       ylogscale=TRUE, fitcol='black', main='!!-plot on danish',
       legendtext=c('gamma', 'Pareto', 'Par-gam', 'Burr'), fitpch=1:4)


##Goodness of Fit Statistics and Tests


data(tplclaimnumber)


x=tplclaimnumber$claim.number

fpois=fitdist(x,'pois')
fnbinom=fitdist(x,'nbinom')

fpoisZM=fitdist(x,'poisZM', start=list(
  prob=sum(x==0)/length(x), lambda=mean(x)),
  lower=c(0,0), upper=c(1,Inf))

gofstat(list(fpois, fnbinom, fpoisZM), chisqbreaks=c(0:4,9),
        discrete=TRUE, fitnames=c('Poisson', 'NegBinomial', 'ZM-Poisson'))


##Skewness-Kurtosis Graph




p=c(.9, .95, .974,.99)
rbind(
  empirical=quantile(danishuni$Loss, prob=p),
  gamma=quantile(fgam, prob=p)$quantiles,
  Pareto=quantile(fpar, prob=p)$quantiles,
  Pareto_gamma=quantile(fmixgampar, prob=p)$quantiles,
  Burr=quantile(fburr,prob=p)$quantiles
)



compmom=function(order)
  c(empirical=sum(danishuni$Loss^order)/length(x),
    gamma=mgamma(order, fgam[[1]][1], fgam[[1]][2]),
    Pareto=mpareto(order, fpar[[1]][1], fpar[[1]][2]),
    Pareto_gamma=as.numeric(fmixgampar[[1]][2]),
    mgamma(order,fmixgampar[[1]][2], fmixgampar[[1]][3])+
      (1-fmixgampar[[1]][1]*
         mpareto(order, fmixgampar[[1]][4], fmixgampar$estimate[5])),
    Burr=mburr(order, fburr[[1]][1],fburr[[1]][2],fburr[[1]][3]))

rbind(Mean=compmom(1), Mom2nd=compmom(2))


descdist(danishuni$Loss, boot=100)

# descdist(tplclaimnumber$claim.number,method='unbiased', discrete = TRUE, boot=500)

var(tplclaimnumber$claim.number)
##이거 왜이렇게 나오는지 모르겠다 


# descdist(data, discrete = FALSE, boot = NULL, method = "unbiased",
#          graph = TRUE, obs.col = "darkblue", obs.pch = 16, boot.col = "orange")
head(tplclaimnumber)



##Linear Regression:
## Introducing Covariates in statistical Inference

#Using Covaiates in the statistical Framework

# library(CASdatasets)
data("Davis")
X=Davis$height

(param=fitdistr(X,'normal')$estimate)

logdf=function(x,parameter){
  p=parameter[1]
  m1=parameter[2]
  s1=parameter[4]
  m2=parameter[3]
  s2=parameter[5]
  return(log(p*dnorm(x,m1,s1)+(1-p)*dnorm(x,m2,s2)))
}

logL=function(parameter) -sum(logdf(X,parameter))

Amat=matrix(c(1,-1,0,0,0,0,
              0,0,0,0,1,0,0,0,0,0,0,0,0,1),4,5)
bvec=c(0,-1,0,0)
constrOptim(c(.5,160,180,10,10), logL, NULL, ui=Amat, ci=bvec)$par

# install.packages('mixtools')

library(mixtools)
mix=normalmixEM(X)
(param12=c(mix$lambda[1], mix$mu,mix$sigma))

sex= Davis$sex
(pM= mean(sex=='M'))

(paramF= fitdistr(X[sex=='F'], 'normal')$estimate)

(paramM= fitdistr(X[sex=='M'], 'normal')$estimate)


f1=function(x) dnorm(x, param[1], param[2])
f2=function(x) {param12[1]*dnorm(x,param12[2], param12[4])+
    (1-param12[1])*dnorm(x,param12[3], param12[5])}
f3=function(x){ pM*dnorm(x,paramM[1], paramM[2]+
                           (1-pM)*dnorm(x,paramF[1], paramF[2]))}
boxplot(X~sex,horizontal=TRUE,names=c('Female', 'Male'))
x=seq(min(X), max(X), by=.1)
plot(x,f2(x), col='grey', lwd=2)
lines(x,f3(x), col='black', lwd=2)
lines(density(X))

#Linear Regression Model

Y=Davis$height
X1=Davis$sex
X2=Davis$weight
df=data.frame(Y,X1,X2)

#Inference in a Linear Model


lin.mod=lm(Y~X1+X2,data=df)
summary(lin.mod)

new.obs=data.frame(X1=c('M', 'M', 'F'), X2=c(100,70,65))
predict(lin.mod, newdata=new.obs)

##Aggregate Loss Distribution


#Computation of the Aggregate Loss Distribution

pgamsum=function(x,dfreq, argfreq, shape, rate, Nmax=10)
{
  tol=1e-10;maxit=10
  nbclaim=0:Nmax
  dnbclaim=do.call(dfreq, c(list(x=nbclaim), argfreq))
  psumfornbclaim=sapply(nbclaim, function(n)
    pgamma(x,shape=shape*n, rate=rate))
  psumtot=psumfornbclaim %*% dnbclaim
  dnbclaimtot=dnbclaim
  iter=0
  while(abs(sum(dnbclaimtot)-1)>tol&&iter<maxit)
  {
    nbclaim=nbclaim+Nmax
    dnbclaim=do.call(dfreq, c(list(x=nbclaim), argfreq))
    psumfornbclaim=sapply(nbclaim, function(n)
      pgamma(x,shape=shape*n,rate=rate))
    psumtot=psumtot+psumfornbclaim%*%dnbclaim
    dnbclaimtot=c(dnbclaimtot, dnbclaim)
    iter=iter+1
  }
  as.numeric(psumtot)
}

parsev=c(3,2); parfreq=10
meansev=mgamma(1,parsev[1],parsev[2])
varsev=mgamma(2,parsev[1],parsev[2])-meansev^2
skewsev=(mgamma(3,parsev[1], parsev[2])-
           3*meansev*varsev-meansev^3)/varsev^(3/2)
meanfreq=varfreq=parfreq[1];skewfreq=1/sqrt(parfreq[1])
meanagg=meanfreq*meansev
varagg=varfreq*(varsev+meansev^2)
skewagg=(skewfreq*varfreq^(3/2)*meansev^3+3*varfreq*meansev*
           varsev+meanfreq*skewsev*varsev^(3/2))/varagg^(3/2)
Fs.s=aggregateDist('simulation',
                   model.freq=expression(y=rpois(parfreq)), model.sev=expression(
                     y=rgamma(parsev[1],parsev[2])),nb.simul=1000)
Fs.n=aggregateDist('normal', moments=c(meanagg,varagg))
Fs.np=aggregateDist('npower', moments=c(meanagg,varagg,skewagg))
Fs.exact=function(x) pgamsum(x,dpois, list(lambda=parfreq),
                             parsev[1],parsev[2], Nmax=100)
x=seq(25,40,length=101)


plot(x,Fs.exact(x), type='l',
     main='Agg. Claim Amount Distribution', ylab='F_s(x)')
lines(x,Fs.s(x), lty=2)
lines(x,Fs.n(x),lty=3)
lines(x,Fs.np(x),lty=4)
legend('bottomright', leg=c('exact', 'simulation',
                            'normal approx.', 'NP approx.'), col='black',
       lty=1:4, text.col='black')



parsev=c(3.1,2*2.1);parfreq=10
xmax=qpareto(1-1e-9, parsev[1], parsev[2])
fx2=discretize(ppareto(x,parsev[1], parsev[2]), from=0,
               to=xmax, step=0.5, method='unbiased',
               lev=levpareto(x,parsev[1],parsev[2]))
Fs2=aggregateDist('recursive',model.freq='poisson',
                  model.sev=fx2, lambda=parfreq, x.scale=0.5, maxit=2000)
fx.u2=discretize(ppareto(x,parsev[1], parsev[2]), from=0,
                 to=xmax, step=0.5, method='upper')
Fs.u2=aggregateDist('recursive', model.freq='poisson',
                    model.sev=fx.u2,lambda=parfreq, x.scale=0.5, maxit=2000)
fx.l2=discretize(ppareto(x,parsev[1], parsev[2]), from=0,
                 to=xmax, step=0.5, method='lower')
Fs.l2=aggregateDist('recursive', model.freq='poisson',
                    model.sev=fx.l2, lambda=parfreq, x.scale=0.5, maxit=2000)





# plot(x,Fs2(x), type='l',
#      main='Agg. Claim Amount Distribution', ylab='F_s(x)')
# 
# lines(x,Fs.u2(x),lty=3)
# lines(x,Fs.l2(x),lty=4)
# legend('bottomright', leg=c('exact', 'simulation',
#                             'normal approx.', 'NP approx.'), col='black',
#        lty=1:4, text.col='black')

#다른 예시는코드가 없다 복원이 필요




##Poisson Process


rate=1
rFexp=function() rexp(1,rate)
rRenewal=function(Tmax=1, rF=rFexp){
  t=0
  vect.W=NULL
  while(t<Tmax){
    W=rF()
    t=t+W
    if(t<T) vect.W=c(vect.W,W)
  }
  
  return(list(T=cumsum(vect.W),W=vect.W,N=length(vect.W)))
}

set.seed(1)
rRenewal(Tmax=2)


rPoissonProc=function(Tmax=1, lambda=rate){
  N=rpois(n=1,lambda*Tmax)
  vect.T=NULL
  if(N>0) vect.T=sort(runif(N))*lambda*Tmax
  return(list(T=vect.T,W=diff(c(0,vect.T)),N=N))
}

set.seed(1)
rPoissonProc(T=5)


lambda=function(t) 100*(sin(t*pi)+1)

Lambda=function(t) integrate(f=lambda, lower=0, upper=t)$value


# install.packages('Ryacas')
library(Ryacas)

??Ryacas
# 
# Tmax=3*pi
# set.seed(1)
# t=1e-15;x=1e-15
# while(X[length(X)]<= Tmax){
#   Ft=function(x) 1-exp(-Lambda(t+x)+Lambda(t))
#   x=uniroot(function(x) Ft(X)-runif(1), interval=c(0,Tmax))$root
#   t=t+x
#   X=c(X,t)}
# X=X[-which.max(X)]
# X
# 극한값을 이용하는거같은데 0에러 발생
# ryacas 라이브러리를 활용할수 있음 시도해보자

lambda.up=200
set.seed(1)
t=0;X=t
while(X[length(X)]<=Tmax){
  u=runif(1)
  t=t-log(u)/lambda.up
  if(runif(1)<=lambda(t)/lambda.up)X=c(X,t)
}
X=X[-c(1,which.max(X))]
X

hist(X, breaks=seq(0, 3*pi,by=pi/32), col='grey',
     border='White', xlab='', main='')
lines(seq(0,3*pi,by=.02), lambda(seq(0, 3*pi,by=.02))*pi/32,lwd=2)

## From Poisson Processes to Levy processes

randX=function(n) rexp(n,1)

rComPoissonProc=function(Tmax=1, lambda=rate, rand){
  N=rpois(n=1, lambda*Tmax)
  X=randX(N)
  vect.T=NULL
  if(N>0) vect.T=sort(runif(N))*lambda*T
  return(list(T=vect.T, W=diff(c(0,vect.T)), X=X,N=N))}

set.seed(1)
rComPoissonProc(Tmax=5,rand=randX)


set.seed(1)
compois=rComPoissonProc(Tmax=5, rand=randX)
St=function(t){sum(compois$X[compois$T<t])}

time=seq(0,5,length=501)
plot(time, Vectorize(St)(time), type='s')
abline(v=compois$T,lty=2,col='grey')


n=1000
h=Tmax/n
set.seed(1)
B=c(0,cumsum(rnorm(n,sd=sqrt(h))))

Bt=function(t){B[trunc(n*t/Tmax)+1]}
time=seq(0,5,length=501)
plot(time, Vectorize(Bt)(time),type='s')
plot(time, Vectorize(L)(time))




mu=lambda(rate)
mu
L=function(t) -mu*t+St(t+Bt(t))


##Ruin models

#Asmussen-Rolski

psi=ruin(claims='e', par.claims=list(rate=1/0.6),
         wait='e', par.wait=list(rate=1/0.6616858))

p=c(0.5614, 0.4386)
r=matrix(c(-8.64, 0.101, 1.997, -1.095), 2,2)
lambda=1/(1.1*mphtype(1,p,r))
psi2=ruin(claims='p', par.claims=list(prob=p,rates=r),
          wait='e', par.wait=list(rate=lambda ))


a=(0.4/5+0.6)*lambda
psi3=ruin(claims='p', par.claims=list(prob=p,rates=r),
          wait='e', par.wait=list(rate=c(5*a,a), weights=
                                    c(0.4,0.6)), maxit=225)


plot(psi, from=0, to=50)
plot(psi2, add=TRUE, lty=2)
plot(psi3, add=TRUE, lty=3)
legend('topright', leg=c('EXp-Exp', 'PH-Exp', 'PH-MixExp'), lty=1:3, col='black')

##Copulas and Multivariate Distributions


##Definition of Copulas

##ArChimedean Copulas
##Elliptical Copulas

##application and copula selection
#the Loss ALAE
library(CASdatasets)
data(lossalae)


par(mfrow=c(1,2))
plot(lossalae, log='xy', main='Scatterplot of loss-ALAE')
plot(apply(lossalae,2,rank)/NROW(lossalae),
     main='rank transform of loss-ALAE')

# install.packages('fCopulae')

library(fCopulae)
dnormcop=function(U,param){
  as.numeric(dellipticalCopula(U,rho=param[1], type='norm'))
}
dtcop=function(U, param){
  as.numeric(dellipticalCopula(U,rho=param[1], type='t',
                               param=param[2]))
}
dgumcop=function(U,param){
  as.numeric(devCopula(U,type='gumbel', param=param[1]))
}
dHRcop=function(U,param){
  as.numeric(devCopula(U,type='husler.reiss', param=param[1]))
}
dfrankcop=function(U,param){
  as.numeric(darchmCopula(U,type='5'), alpha=param[1])
}


paretochart=function(x){
  plot(-log((1:length(x))/length(x)+1), log(sort(x)))}

paretochart(lossalae$Loss)
paretochart(lossalae$ALAE)

fit.cop.IFM.2=function(obs, copula, marg, arg.margin=list(),
                       method.margin='mle', arg.cop=list(), initpar){
  Obs1=obs[,1]
  Obs2=obs[,2]
  if(marg %in% c('exp', 'gamma', 'lnorm', 'pareto', 'burr')){
    Obs1=Obs1[Obs1>0]
    Obs2=Obs2[Obs2>0]}
  marg1=do.call(fitdist, c(list(data=Obs1, distr=marg,
                                method=method.margin), arg.margin))
  marg2=do.call(fitdist, c(list(data=Obs2, distr=marg,
                                method=method.margin), arg.margin))
  comput.cdf=function(fit,obs){
    para=c(as.list(fit$estimate), as.list(fit$fix.arg))
    distname=fit$distname
    pdistname=paste('p', distname, sep='')
    do.call(pdistname, c(list(q=obs), as.list(para)))  }
  pseudomarg1=comput.cdf(marg1, Obs1)
  pseudomarg2=comput.cdf(marg2, Obs2)
  U=cbind(pseudomarg1, pseudomarg2)
  copLogL=function(x){
    if(arg.cop$lower<=x && arg.cop$upper>=x)
      res=-sum(remove.naninf(log(copula(U,param=x))))
    else res=Inf
    return(res)  }
  resopt=optim(par=initpar, fn=copLogL, method='L-BFGS-B',
               lower=arg.cop$lower, upper=arg.cop$upper)
  list(marg1=marg1, marg2=marg2, copula=
         list(name=arg.cop$name, alpha=resopt$par))  }
remove.naninf <- function(x) {
  x <- x[which(!is.nan(x) & is.finite(x))]
  return(x)
}




library(fCopulae)
argnorm=list(length=1, lower=0, upper=1, name='Gaussian')
argt=list(length=2, lower=c(0,0), upper=c(1,1000),
          name='Student')
arggum=list(length=1, lower=1, upper=100, name='Gumbel')
argHR=list(length=1, lower=0, upper=1000, name='Husler-Reiss')
argfrank=list(length=1, lower=0, upper=1000, name='Frank')
fgausspareto=fit.cop.IFM.2(lossalae, copula=dnormcop,
                           marg='pareto', arg.margin=list(start=list(shape=10, scale=100),
                                                          lower=c(1,1/2)), arg.cop=argnorm, initpar=1/2)
ftpareto=fit.cop.IFM.2(lossalae, copula=dtcop,
                       marg='pareto', arg.margin=list(start=list(shape=10, scale=100),
                                                      lower=c(1,1/2)), arg.cop=argt, initpar=c(1/2,4))
fgumbelpareto=fit.cop.IFM.2(lossalae, copula=dgumcop,
                            marg='pareto', arg.margin=list(start=list(shape=10, scale=100),
                                                           lower=c(1,1/2)), arg.cop=arggum, initpar=10)
fHRpareto=fit.cop.IFM.2(lossalae, copula=dHRcop,
                        marg='pareto', arg.margin=list(start=list(shape=10, scale=100),
                                                       lower=c(1,1/2)), arg.cop=argHR, initpar=10)
ffrankpareto=fit.cop.IFM.2(lossalae, copula=dfrankcop,
                           marg='pareto', arg.margin=list(start=list(shape=10, scale=100),
                                                          lower=c(1,1/2)), arg.cop=argfrank, initpar=10)


recap=function(x){
  res=c(alpha=x$copula$alpha, x$marg1$estimate, x$marg2$estimate)
  if(length(res)<6)
    res=c(res[1], NA, res[2:5])
  res=as.matrix(res)
  colnames(res)=x$copula$name
  res}

round(cbind(recap(fgausspareto), recap(ftpareto),
            recap(fHRpareto), recap(fgumbelpareto),
            recap(ffrankpareto)),4)                                                          
                              
                                                          

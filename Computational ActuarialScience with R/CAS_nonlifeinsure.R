## non-life insurance

## 일반 보험에서 집단 모형

##표준 보험료 원리는 기대값원리

##Pure Premium in a Heterogenous Context

library(CASdatasets)
data("freMTPLfreq")
contracts=freMTPLfreq
contracts_f=freMTPLfreq
data('freMTPLsev')

# claims=freMTPLsev

contracts_f$DriverAge=cut(contracts$DriverAge,c(17,22,26,42,74,Inf))
contracts_f$CarAge=cut(contracts$CarAge,c(0,1,4,15,Inf), include.lowest=TRUE)
contracts_f$Density=cut(contracts$Density,c(0,40,200,500,4500,Inf),include.lowest=TRUE)
#데이터 로드후 분할점 적용





##Claims Frequency and Log -poisson Regression


##청구건수와 가입기간으로 평균과 분산 파이계산
#파이값은 포아송분포의 파이
vY=contracts_f$ClaimNb
vE=contracts_f$Exposure
m=sum(vY)/sum(vE)
v=sum((vY-m*vE)^2)/sum(vE)
cat("average=",m,"variance=",v,"phi=",v/m,"\n")
# average= 0.06979859 variance= 0.07396742 phi= 1.059727 
#계산결과 파이값은 1에 근접하다
##포아송분포의 특징: 평균과 분산이 같다(이론상)
##따라서 청구건수와 가입기간은 포아송 분포를 따른다


## 포아송분포에서는 파이가 l과 1이어야 한다
#범주별 계산
vX=as.factor(contracts_f$Region)
for(i in 1:length(levels(vX))){
  vEi=vE[vX==levels(vX)[i]]
  vYi=vY[vX==levels(vX)[i]]
  mi=sum(vYi)/sum(vEi)
  vi=sum((vYi-mi*vEi)^2)/sum(vEi)
  cat('average=', mi, 'variance=', vi, 'phi=',vi/mi,'\n')
  
}
## 지역과 상관없이 파이가 1
# average= 0.07366204 variance= 0.08064593 phi= 1.09481 
# average= 0.06788926 variance= 0.07152207 phi= 1.053511 
# average= 0.06741501 variance= 0.07016521 phi= 1.040795 
# average= 0.06303979 variance= 0.06483763 phi= 1.028519 
# average= 0.06923285 variance= 0.07388641 phi= 1.067216 
# average= 0.08577016 variance= 0.09526731 phi= 1.110728 
# average= 0.08222055 variance= 0.08952784 phi= 1.088874 
# average= 0.08210142 variance= 0.09201134 phi= 1.120703 
# average= 0.07185182 variance= 0.07590899 phi= 1.056466 
# average= 0.0716627 variance= 0.07559456 phi= 1.054866 
# 

## 포아송 회귀
##최대 우도법사용
Y=contracts$ClaimNb
E=contracts$Exposure
(lambda=sum(Y)/sum(E))
#[1] 0.06979859
#람다(평균)계산
weighted.mean(Y/E,E)
#가중평균 계산
# [1] 0.06979859
dpois(0:3,lambda)*100
# [1] 93.258163240  6.509288286  0.227169572  0.005285372
#dpois(x, lambda, log = FALSE)
# 설명-평균(=분산)이 람다인 포아송분포의 확률밀도함수
# x를 0,1,2,3순차적으로 대입하여 확률값 반환
# 예 사건수 0이하로 발생할 확률= dpois(0,lambda)*100

#glm(y~X1+X2+X3+offset(E), family=poisson(link='log'))
#로그 포아송 glm

reg=glm(ClaimNb~Gas+DriverAge+Density+offset(log(Exposure)), family=poisson,
        data=contracts_f)


# summary(reg)
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.6868  -0.3469  -0.2687  -0.1511   6.4723  
#잔차에 대한 분위수 

# Coefficients: 계수

#   변수            표준편차      표준편차  z값       p값
#                Estimate     Std. Error  z value     Pr(>|z|)    
# (Intercept)   -2.144e+00    2.706e-02 -79.244   <2e-16 ***
#   GasRegular  -1.350e-01  1.592e-02  -8.479   <2e-16 ***
#   DriverAge   -1.069e-02  5.623e-04 -19.014   <2e-16 ***
#   Density      2.179e-05  1.478e-06  14.740   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 105613  on 413168  degrees of freedom
# Residual deviance: 104969  on 413165  degrees of freedom
# AIC: 136234
# # 
# # Number of Fisher Scoring iterations: 6
# # This   stucture   of   the   output   is   very   close   to 
# the   one   obtained   with   the   linear   regression 
# # function  lm,  in  R.  See  Chambers  &  Hastie  (1991)  and  Faraway  (2006)
# for  more  details.
# # In   the   following   sections,   we   will   discuss   the   interpretation 
# of   this   regression,   on   categorical  variables  (one  or  two)  and 
# on  continuous  variables  (one  or  two).

# Ratemaking   with   One   Categorical   Variable

#단일 범주형 변수 ratemaking
#Gas type으로 분석
vY=contracts_f$ClaimNb
vE=contracts_f$Exposure
X1=contracts_f$Gas
name1=levels(X1)
##연간 청구횟수


tapply(vY,X1,sum)
# Diesel Regular 
# 8446    7735
# Gastype 별 합계

tapply(vY,X1,sum)/tapply(vE, X1, sum)
# Diesel    Regular 
# 0.07467412 0.06515364
# Gastype 별 연간 청구횟수


df=data.frame(vY,vE,X1)
regpoislog=glm(vY~0+X1+offset(log(vE)),data=df, family=poisson(link='log'))
summary(regpoislog)

##로그화된 포아송회귀 
# glm(formula = vY ~ 0 + X1 + offset(log(vE)), family = poisson(link = "log"), 
#     data = df)

# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.5092  -0.3610  -0.2653  -0.1488   6.5858  
# 
# Coefficients:
#           Estimate Std. Error z value Pr(>|z|)    
# X1Diesel  -2.59462    0.01088  -238.5   <2e-16 ***
#   X1Regular -2.73101    0.01137  -240.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 450747  on 413169  degrees of freedom
# Residual deviance: 105537  on 413167  degrees of freedom
# AIC: 136799
# 
# Number of Fisher Scoring iterations: 6
# 회귀 모형 적합하다

exp(coefficients(regpoislog))


newdf=data.frame(X1=name1, vE=rep(1,length(name1)))
#예측기능 수행

predict(regpoislog, newdata=newdf, type="response")

regpoislog=glm(vY~X1+offset(log(vE)), data=df,
               family=poisson(link='log'))
summary(regpoislog)



# Call:
#   glm(formula = vY ~ X1 + offset(log(vE)), family = poisson(link = "log"), 
#       data = df)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.5092  -0.3610  -0.2653  -0.1488   6.5858  
# 
# Coefficients:
#               Estimate Std. Error z value    Pr(>|z|)    
# (Intercept) -2.59462    0.01088 -238.454   <2e-16 ***
#   X1Regular   -0.13639    0.01574   -8.666   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 105613  on 413168  degrees of freedom
# Residual deviance: 105537  on 413167  degrees of freedom
# AIC: 136799
# 
# Number of Fisher Scoring iterations: 6
#모형적합 확인

exp(coefficients(regpoislog))
# (Intercept)   X1Regular 
# 0.07467412  0.87250624 
## 잔차와 x1  총 청구건수중 디젤차는 0.0746 -> 87.25%는 휘발유차?


prod(exp(coefficients(regpoislog)))






exp(coefficients(regpoislog))
# X1Diesel  X1Regular 
# 0.07467412 0.06515364 
# exponential of the coefficents (로그화된) 회귀계수 를 지수화

prod(exp(coefficients(regpoislog)))
# [1] 0.004865291
# 로그화된 회귀계수 곱
#왜 곱을 구했나?
#


# Contingency   Tables   and   Minimal   Bias   Techniques
# 분할표

X2=contracts_f$Density
names1=levels(X1)
names2=levels(X2)
(P=table(X1,X2))
## 연료타입과 인구밀도(운전자기준)로 분할표 작성

E=Y=P

for(k in 1: length(names1)){
  E[k,]=tapply(vE[X1==names1[k]],X2[X1==names1[k]],sum)
  Y[k,]=tapply(vY[X1==names1[k]], X2[X1==names1[k]], sum)
}
#익스포져 매트릭스 E, 청구건수 매트릭스 Y생성

E

Y

(N=Y/E)

L=matrix(NA,100,length(names1))
C=matrix(NA,100,length(names2))

C[1,]=rep(sum(vY)/sum(vE),length(names2));colnames(C)=names2

for(j in 2:100){
  for(k in 1:length(names1)) L[j,k]=sum(Y[k,])/sum(E[k,]*C[j-1,])
  for(k in 1:length(names2)) C[j,k]=sum(Y[,k])/sum(E[,k]*L[j,])
}


L[100,]
C[100,]

PredN=N
for(k in 1:length(names1)) PredN[k,]=L[100,k]*C[100,]
PredN



sum(PredN[1,]*E[1,])


sum(Y[1,])


df=data.frame(vY,vE,X1,X2)
regpoislog=glm(vY~X1+X2, offset=log(vE), data=df,
family=poisson(link='log'))
newdf=data.frame(
  X1=factor(rep(names1,length(names2))),
  vE=rep(1,length(names1)*length(names2)),
  X2=factor(rep(names2,each=length(names1))))


matrix(predict(regpoislog,newdata=newdf,
               type='response'),length(names1),length(names2))

reg.cut=glm(ClaimNb~DriverAge+offset(log(Exposure)), family=poisson,
            data=contracts_f)

summary(reg.cut)

reg.poisson=glm(ClaimNb~DriverAge+offset(log(Exposure)), family=poisson,
                data=contracts)
summary(reg.poisson)

newdb=data.frame(DriverAge=18:99,Exposure=1)
pred.poisson=predict(reg.poisson,newdata=newdb, type='response', se=TRUE)
plot(18:99,pred.poisson$fit, type='l', xlab='Age of the driver',
     ylab='Annualized Frequency', ylim=c(0,.3), col='grey',lwd=7)
segments(18:99,pred.poisson$fit-2*pred.poisson$se.fit,
         18:99,pred.poisson$fit+2*pred.poisson$se.fit,col='grey',lwd=7)
lines(18:99,pred.poisson$fit)
abline(h=sum(contracts$ClaimNb)/sum(contracts$Exposure),lty=2)

library(actuar)

discretize()

reg.np=glm(ClaimNb~as.factor(DriverAge)+offset(log(Exposure)), family=poisson,
                data=contracts)
summary(pred.poissonnp)

pred.poissonnp=predict(reg.np,newdata=newdb, type='response', se=TRUE)
# discpred=discretize(poisson()  pred.poissonnp)

plot(18:99,pred.poisson$fit, type='l', xlab='Age of the driver',
     ylab='Annualized Frequency', ylim=c(0,.3), col='grey',lwd=7)
segments(18:99,pred.poissonnp$fit-2*pred.poissonnp$se.fit,
         18:99,pred.poissonnp$fit+2*pred.poissonnp$se.fit,col='grey',lwd=7)
lines(18:99,pred.poissonnp$fit)
abline(h=sum(contracts$ClaimNb)/sum(contracts$Exposure),lty=2)
#이산화는 일단 넘기자... 시간이 없다


library(nlme)
library(mgcv)
reg.splines=gam(ClaimNb~s(DriverAge)+offset(log(Exposure)),
                family=poisson, data=contracts)

newdb=data.frame(DriverAge=18:99,Exposure=1)
pred.reg=predict(reg.splines,newdata=newdb, type='response', se=TRUE)
plot(18:99,pred.reg$fit, type='l', xlab='Age of the driver',
     ylab='Annualized Frequency', ylim=c(0,.3), col='grey',lwd=7)
segments(18:99,pred.reg$fit-2*pred.reg$se.fit,
         18:99,pred.reg$fit+2*pred.reg$se.fit,col='grey',lwd=7)
lines(18:99,pred.reg$fit)
abline(h=sum(contracts$ClaimNb)/sum(contracts$Exposure),lty=2)
##일반화 가법모형
summary(reg.splines)


contracts_f$CarAge=cut(contracts$CarAge,c(0,15,Inf), include.lowest=TRUE)
levels(contracts_f$CarAge)

contracts_f$Powerf=factor(1*(contracts_f$Power%in%letters[4:6])+
                            2*(contracts_f$Power%in%letters[7:8]),
                          labels=c('other','DEF','GH'))



contracts_f$Brandf=factor(contracts_f$Brand!='other',labels=c('other','F'))
factor(contracts_f$Brandf)



freq=formula(ClaimNb~DriverAge+CarAge+Density+Brandf+Powerf+Gas+offset(log(Exposure)))
regp=glm(freq, data=contracts_f, family=poisson(link='log'))
summary(regp)


library(MASS)
regnb2=glm.nb(freq, data=contracts_f)
summary(regnb2)
##콰시 포아송  음이항분포1형 2형이있는데 1형은 편향이 심하다 그래서 2형으로

# 
# library(MASS)
# glm.nb(Y~X1+X2+X3+offset(log(E)))
# glm(Y~X1+X2+X3+offset(log(E)), family=negative.binomial(1))
#이거 왜 설명했는지 모르겠다

install.packages('AER')


library(AER)
dispersiontest(regp)


# zero inflatedmodel
# 영과잉 모형
install.packages('gamlss')
library(splines)
library(gamlss.data)
library(gamlss)
library(gamlss.zip)

fregzi=formula(ClaimNb~DriverAge+CarAge+Density+Brandf+Powerf+
                 Gas+offset(log(Exposure))|1)

regzip=zeroinfl(fregzi, data=contracts_f, dist='poisson',
                link='logit')

summary(regzip)



fregzi=formula(ClaimNb~DriverAge+CarAge+Density+Brandf+Powerf+
                 Gas+offset(log(Exposure))|DriverAge)

regzip=zeroinfl(fregzi, data=contracts_f, dist='poisson', link='logit')

claims=freMTPLsev
claims_f=freMTPLsev

claims=merge(claims,contracts)
claims_f=merge(claims_f,contracts_f)





reg.logn=lm(log(ClaimAmount)~CarAge+Gas,
            data=claims_f[claims_f$ClaimAmount<15000,])

reg.gamma=glm(ClaimAmount~CarAge+Gas, family=Gamma(link='log'),
              data=claims_f[claims_f$ClaimAmount<15000,])

summary(reg.gamma)
summary(reg.logn)




reg.logn=lm(log(ClaimAmount)~DriverAge, data=claims)
reg.gamma=glm(ClaimAmount~DriverAge, family=Gamma(link='log'), data=calims)
summary(reg.gamma)



#Large Claims and Ratemaking

reg.logn=lm(log(ClaimAmount)~DriverAge,data=claims)
reg.gamma=glm(ClaimAmount~DriverAge, family=Gamma(link='log'), data=claims)
summary(reg.gamma)


summary(reg.logn)

mean(claims$ClaimAmount)

mean(predict(reg.gamma,type='response'))

sigma=summary(reg.logn)$sigma
mean(exp(predict(reg.logn))*exp(sigma^2/2))

M=claims[order(-claims$ClaimAmount),c('ClaimAmount','ClaimNb','Power',
         'CarAge','DriverAge', 'Gas', 'Density')]
M$sum=cumsum(M$ClaimAmount)/sum(M$ClaimAmount)

head(M)

s=10000

claims$Standard=(claims$ClaimAmount<s)

mean(claims$Standard)


library(splines)
age=seq(18,100)
regC=glm(Standard~bs(DriverAge), data=claims, family=binomial)
ypC=predict(regC, newdata=data.frame(DriverAge=age),type='response',
            se=TRUE)


plot(age,ypC$fit,ylim=c(.95,1),type='l',)
polygon(c(age,rev(age)),c(ypC$fit+2*ypC$se.fit,rev(ypC$fit-2*ypC$se.fit)),
        col=adjustcolor('grey', alpha=0.2),border=NA)
abline(h=mean(claims$Standard),lty=2)

indexstandard=which(claims$ClaimAmount<s)

mean(claims$ClaimAmount[indexstandard])
mean(claims$ClaimAmount[-indexstandard])


regA=glm(ClaimAmount~bs(DriverAge),data=claims[indexstandard,],
         family=Gamma(link='log'))
ypA=predict(regA,newdata=data.frame(DriverAge=age),type='response')

regB=glm(ClaimAmount~bs(DriverAge), data=claims[-indexstandard,],
         family=Gamma(link='log'))
ypB=predict(regB,newdata=data.frame(DriverAge=age), type='response')


reg=glm(ClaimAmount~bs(DriverAge), data=claims,family=Gamma(link='log'))
yp=predict(reg,newdata=data.frame(DriverAge=age), type='response')


ypC=predict(regC,newdata=data.frame(DriverAge=age), type='response')

plot(age,yp,type='l',lwd=2, ylab='Average cost',xlab='Age of the driver')
lines(age,ypC*ypA+(1-ypC)*ypB,type='h',col='grey', lwd=6)
lines(age,ypC*ypA, type='h', col='black',lwd=6)
abline(h=mean(claims$ClaimAmount), lty=2)



install.packages('nnet')
library(mgcv)
library(nnet)

threshold=c(0,1150,10000,Inf)
regD=multinom(cut(claims$ClaimAmount, breaks=threshold)~bs(DriverAge),data=claims)

summary(regD)

A=tapply(claims$ClaimAmount, claims$PolicyID,sum )

ADF=data.frame(PolicyID=names(A),ClaimAmount=as.vector(A))
CT=merge(contracts,ADF,all.x=TRUE)
CT$ClaimAmount[is.na(CT$ClaimAmount)]=0
tail(CT)

CT_f=merge(contracts_f,ADF,all.x=TRUE)


CT_f$ClaimAmount[is.na(CT_f$ClaimAmount)]=0


tail(CT_f)

# install.packages('tweedie')
# install.packages('statmod')
library(tweedie)
library(statmod)
out=tweedie.profile(ClaimAmount~Power+CarAge+
                      DriverAge+Brand+Gas+Density,data=CT_f,
                    p.vec=seq(1.05,1.95,by=.05))
#주의:연산이 오래걸린다


# warnings()

out$p.max

plot(out,type='b')
abline(v=out$p.max,lty=2)

reg1=glm(ClaimAmount~Power+CarAge+DriverAge+Brand+Gas+Density,
         data=CT_f, family=tweedie(var.power = 1, link.power=0),
         start=reg1$coefficients)$coefficients}

#왜 뜬금없이 중괄호가 나오는거지??

vp=seq(1,2,by=.1)
Cp=Vectorize(coef)(vp)
matplot(vp,t(Cp[-1,]),type='1')
text(2,Cp[-1,length(vp)], rownames(Cp[-1,]),cex=.5,pos=4)


install.packages('cplm')


# 여기까지 General insurance pricing 끝




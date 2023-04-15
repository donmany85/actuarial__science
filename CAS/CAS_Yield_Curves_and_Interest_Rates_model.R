library(CASdatasets)

data(FedYieldCurve)

FedYieldCurve

require(xts)

# install.packages('YieldCurve')

require(YieldCurve)

first(FedYieldCurve,'3 month')
last(FedYieldCurve,'3 month')

mat.Fed=c(3/12, 0.5, 1,2,3,5,7,10)
par(mfrow=c(2,3))
for (i in c(1,2,3,370,371,372)){
  plot(mat.Fed, FedYieldCurve[i,], type='o',
       xlab='Maturities structure in years', ylab='Interest rates values')
  title(main=paste('Federal Reserve yield curve observed at',
                   time(FedYieldCurve[i], sep=' ')))
  grid()
}

FedYildCurve$maturity

mat.Fed

seq(mat.Fed)

length(mat.Fed)

head(FedYieldCurve)


par(mfrow=c(1,1))
persp(1982:2012, c(3/12, 0.5, 1,2,3,5,7,10),
      as.matrix(FedYieldCurve[seq(2,nrow(FedYieldCurve),by=12),]),
      theta=30, xlab='year', ylab='Maturity(in years)',
      zlab='Interest rates (in%)', ticktype='detailed', shade=0.2, expand=.3)


M=as.matrix(FedYieldCurve)

(pca.rates=princomp(M,scale=TRUE))
summary(pca.rates)


factor.loadings=pca.rates$loadings[,1:3]
matplot(mat.Fed, factor.loadings, type='l', lwd=c(2,1,1),
        lty=c(1,1,2), xlab='Maturity (in years)', ylab='Factor loadings')

##알아봤는데 맞는 그래프다 eigen을 이용할경우 예시대로 나오긴 한다
##https://www.r-bloggers.com/2022/05/understanding-pca-3-factors-of-the-yield-curve-using-r-code/



##Nelson-Siegel model

factorBeta1=function(lambda, Tau)
{
  (1-exp(-Tau/lambda))/(Tau/lambda)
}
maturity.set=c(3/12,6/12,seq(1:30))
lambda.set=c(0.5,1,2,3,4,5)
par(mfrow=c(2,3))
for( i in 1: length(lambda.set)){
  FB1=factorBeta1(lambda.set[i],maturity.set)
  plot(maturity.set, FB1, type='o', ylim=c(0,1))
  text(12,0.9, substitute(list(lambda)==group('',list(x),''), list(x=i)),cex=1.5)
}




factorBeta2=function(lambda, Tau)
{(1-exp(-Tau/lambda))/(Tau/lambda)-exp(-Tau/lambda)
}
par(mfrow=c(2,3))
for( i in 1: length(lambda.set)){
  FB2=factorBeta2(lambda.set[i], maturity.set)
  plot(maturity.set, FB2, type='o', ylim=c(0,0.4))
  text(i+2,0.35, substitute(list(lambda)==group('',list(x),''),list(x=i)),
       cex=1.5)
  abline(v=i,lty=2)
}



Fed.Rate1=Nelson.Siegel(first(FedYieldCurve,'3 month'),mat.Fed)

Fed.Rate1


Fed.Rate2=Nelson.Siegel(last(FedYieldCurve, '3 month'), mat.Fed)

Fed.Rate2


system.time(Nelson.Siegel(FedYieldCurve, mat.Fed))

Fed.Rates=Nelson.Siegel(FedYieldCurve, mat.Fed)

first(Fed.Rates, 'year')
last(Fed.Rates,'year')


par(mfrow=c(2,2))
plot(Fed.Rates$beta_0, main='Beta_0 coefficient', ylab='Values')

plot(Fed.Rates$beta_1, main='Beta_1 coefficient', ylab='Values')

plot(Fed.Rates$beta_2, main='Beta_2 coefficient', ylab='Values')

plot(Fed.Rates$lambda, main='Lambda coefficient', ylab='Values')


last(Fed.Rates, 'month')



b0=rep(6.590762,21)
b1=rep(-6.49626,21)
b2=seq(-10,10,by=1) #create a sequence of fictive beta_2 coeff
lambda=rep(0.1839190,21)
B=ts(cbind(b0,b1,b2,lambda), start=c(2000,1,1),
     frequency=12) #create a time series object
B=as.xts(B,RECLASS=FALSE) # transform the time series object in xts
A=NSrates(B,mat.Fed) #create the fictive yield curves



#create an interactive plot shows the movement of the yield curve
#at different values of beta_2

for (i in 1:nrow(A))
{
  plot(mat.Fed, A[i,],type='l', ylim=c(-1,6))
  title(main=paste('beta_2', B[i,3], sep='='))
  par(ask=TRUE)
}
  


##Svensson model
## 코드가 없다 2000년대에 소개되었는데 성능은 좋다고 한다
## Nelson-Sigel model 과 유사하다고 한다








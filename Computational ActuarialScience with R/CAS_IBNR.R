##Claims Reserving and IBNR


library(CASdatasets)
data(ukaggclaim)

#chain ladder

n=7
Claims=data.frame(originf=factor(rep(2007:2013, n:1)),
                  dev=sequence(n:1),
                  inc.paid=
                    c(3511, 3215, 2266, 1712, 1059, 587,
                      
                      340, 4001, 3702, 2278, 1180,
                      956,
                      
                      629, 4355, 3932, 1946, 1522, 1238,
                      
                      4295, 3455, 2023, 1320, 4150, 3747,
                      
                      2320, 5102, 4548, 6283))


#make triangle
(inc.triangle=with(Claims, {
  M=matrix(nrow=n, ncol=n,
           dimnames=list(origin=levels(originf), dev=1:n))
  M[cbind(originf, dev)]=inc.paid
  M
}))


#cumalative
(cum.triangle=t(apply(inc.triangle, 1, cumsum)))


(latest.paid=cum.triangle[row(cum.triangle)==n-col(cum.triangle)+1])


Claims$cum.paid=cum.triangle[with(Claims, cbind(originf, dev))]

op=par(fig=c(0,0.5,0,1), cex=0.8, oma=c(0,0,0,0))
with(Claims, {
  interaction.plot(x.factor=dev, trace.factor=originf, response=inc.paid,
                   fun=sum, type='b', bty='n', legend=FALSE);
  axis(1,at=1:n)
  par(fig=c(0.45, 1,0,1), new=TRUE, cex=0.8, oma=c(0,0,0,0))
  interaction.plot(x.factor=dev, trace.factor=originf, response=cum.paid,
                   fun=sum, type='b', bty='n');axis(1, at=1:n)})

mtext('Incremental and cumulative claims development',
      side=3, outer=TRUE, line=-3, cex=1.1, font=2)
par(op)


library(lattice)
xyplot(cum.paid~dev|originf, data=Claims, t='b', layout=c(4,2),
       as.table=TRUE, main='Cumulative claims development')





f=sapply((n-1):1, function(i){
  sum(cum.triangle[1:i, n-i+1])/sum(cum.triangle[1:i, n-i])})

length(f)

tail=1
(f=c(f,tail))


full.triangle=cum.triangle
for( k in 1:(n-1)){
  full.triangle[(n-k+1):n, k+1]=full.triangle[(n-k+1):n,k]*f[k]
}
full.triangle

(ultimate.paid=full.triangle[,n])



(ldf=rev(cumprod(rev(f))))

(dev.pattern=1/ldf)

(reserve=sum(latest.paid*(ldf-1)))
#or
sum(ultimate.paid-latest.paid)


a=ultimate.paid 
(b=c(dev.pattern[1], diff(dev.pattern)))

(X.hat=a%*%t(b))


(BF2013=ultimate.paid[n]*dev.pattern[1]+20000*(1-dev.pattern[1]))



#tail factors

dat=data.frame(lf1=log(f[-c(1,n)]-1), dev=2:(n-1))
(m=lm(lf1~dev, data=dat))


par("mar")

par(mar=c(1,1,1,1))
plot(lf1~dev, main='log(f-1)~dev', data=dat, bty='n')
abline(m)

sigma=summary(m)$sigma
extrapolation=predict(m, data.frame(dev=n:100))
(tail=prod(exp(extrapolation+0.5*sigma^2)+1))

# install.packages('ChainLadder')
library(ChainLadder)

ata(cum.triangle)


#stochastic reserving models 

names(Claims)[3:4]=c('inc.paid.k', 'cum.paid.k')
ids=with(Claims, cbind(originf, dev))
Claims=within(Claims,{
  cum.paid.kp1=cbind(cum.triangle[,-1], NA)[ids]
  inc.paid.kp1=cbind(inc.triangle[,-1], NA)[ids]
  devf=factor(dev)
})


delta=0:2 
ATA=sapply(delta, function(d)
  coef(lm(cum.paid.kp1~0+cum.paid.k:devf, weights=1/cum.paid.k^d, data=Claims)))
dimnames(ATA)[[2]]=paste('Delta=',delta)
ATA


xyplot(cum.paid.kp1~cum.paid.k|devf, 
       data=subset(Claims, dev<(n-1)),
       main='Age-to-age developments', as.table=TRUE,
       scales=list(relation='free'),
       key=list(columns=1,lines=list(lty=1:4, type='l'),
                text=list(lab=c('lm(y~x)',
                                'lm(y~0+x)',
                                'lm(y~0+x,w=1/x)',
                                'lm(y~0+x, w=`/x^2'))),
       panel=function(x,y,...){
         panel.xyplot(x,y,...)
         if(length(x)>1){
           panel.abline(lm(y~x),lty=1)
           panel.abline(lm(y~0+x),lty=2)
           panel.abline(lm(y~0+x, weights=1/x), lty=3)
           panel.abline(lm(y~0+x,,weights=1/x^2),lty=4)
         }
       })


#Mack model

#library(ChinLadder)

(mack=MackChainLadder(cum.triangle, weight=1, alpha=1, est.sigma='Mack'))

# MackChainLadder(Triangle=cum.triangle, weights=1, alpha=1,
#                est.sigma='Mack')


plot(mack, lattice=TRUE, layout=c(4,2))
plot(mack)


#Poisson Regression Model for Incremental Claims 

preg=glm(inc.paid.k~originf+devf,
         data=Claims, family=poisson(link='log'))

summary(preg)

allClaims=data.frame(origin=sort(rep(2007:2013, n)),
                     dev=rep(1:n,n))
allClaims=within(allClaims, {
  devf=factor(dev)
  cal=origin+dev-1 
  originf=factor(origin)
})
(pred.inc.tri=t(matrix(predict(preg,type='response',
                               newdata=allClaims), n,n)))

df=c(0,coef(preg)[(n+1):(2*n-1)])
sapply(2:7, function(i) sum(exp(df[1:i]))/sum(exp(df[1:(i-1)])))


# install.packages('AER')
library(AER)              
dispersiontest(preg)

odpreg=glm(inc.paid.k~originf+devf, data=Claims,
           family=quasipoisson)
summary(odpreg)



mu.hat=predict(odpreg,newdata=allClaims, type='response')*(allClaims$cal>2013)
phi=summary(odpreg)$dispersion 
Sigma=vcov(odpreg)
model.formula=as.formula(paste('~', formula(odpreg)[3]))


#future design matrix 

X=model.matrix(model.formula, data=allClaims)

Cov.eta=X%*%Sigma%*%t(X)

sqrt(phi*sum(mu.hat)+t(mu.hat)%*%Cov.eta%*%mu.hat)       

op=par(mfrow=c(2,2), oma=c(0,0,3,0))
plot(preg)
par(op)


(odp=glmReserve(as.triangle(inc.triangle), var.power=1, cum=FALSE))

#bootstrap chainladder

set.seed(1)
(B=BootChainLadder(cum.triangle, R=1000, process.distr='od.pois'))

plot(B)

quantile(B, c(0.75,0.95,0.99,0.995))


# install.packages('fitdistrplus')
library(fitdistrplus)
(fit=fitdist(B$IBNR.Totals[B$IBNR.Totals>0], 'lnorm'))

plot(fit)

qlnorm(0.995, fit$estimate['meanlog'], fit$estimate['sdlog'])


ny=(col(inc.triangle)==(nrow(inc.triangle)-row(inc.triangle)+2))
paid.ny=apply(B$IBNR.Triangles, 3,
              function(x){
                next.year.paid=x[col(x)==(nrow(x)-row(x)+2)]
                sum(next.year.paid)
              })

paid.ny.995=B$IBNR.Triangles[,,order(paid.ny)[round(B$R*0.995)]]
inc.triangle.ny=inc.triangle
(inc.triangle.ny[ny]=paid.ny.995[ny])


Claims=within(Claims, {
  log.inc=log(inc.paid.k)
  cal=as.numeric(levels(originf))[originf]+dev-1
})

Claims=within(Claims, {
  d1=ifelse(dev<2,1,0)
  d27=ifelse(dev<2,0,dev-1)
})
fit1=lm(log.inc~originf+d1+d27,data=Claims)
summary(fit1)


Claims=within(Claims, {
  a6=ifelse(originf==2012,1,0)
  a7=ifelse(originf==2013,1,0)
  })
fit2=lm(log.inc~a6+a7+d1+d27,data=Claims)
summary(fit2)


op=par(mfrow=c(2,2),oma=c(0,0,3,0))
plot(fit2)
par(op)


shapiro.test(fit2$residuals)

resPlot=function(model,data){
  xvals=list(
    fitted=model[['fitted.values']],
    origin=as.numeric(levels(data$originf))[data$originf],
    cal=data$cal, dev=data$dev
  )
op=par(mfrow=c(2,2), oma=c(0,0,3,0))
for (i in 1:4){
  plot.default(rstandard(model)~xvals[[i]],
               main=paste('Residuals vs', names(xvals)[i]),
               xlab=names(xvals)[i], ylab='Standardized residuals')
  panel.smooth(y=rstandard(model), x=xvals[[i]])
  abline(h=0, lty=2)
}
mtext(as.character(model$call)[2], outer=TRUE, cex=1.2)
par(op)
}

resPlot(fit2,Claims)


Claims=within(Claims, {
  p34=ifelse(cal<2011&cal>2008, cal-2008,0)
})
fit3=update(fit2,.~ .+p34,data=Claims)
summary(fit3)

resPlot(fit3,Claims)


log.incr.predict=function(model, newdata){
  Pred=predict(model, newdata=newdata, se.fit=TRUE)
  Y=Pred$fit 
  VarY=Pred$se.fit^2+Pred$residual.scale^2
  P=exp(Y+VarY/2)
  VarP=P^2*(exp(VarY)-1)
  seP=sqrt(VarP)
  model.formula=as.formula(paste('~', formula(model)[3]))
  mframe=model.frame(model.formula, data=newdata)
  X=model.matrix(model.formula, data=newdata)
  varcovar=X%*%vcov(model) %*%t(X)
  CoVar=sweep(sweep((exp(varcovar)-1),1,P,'*'),2,P,'*')
  CoVar[col(CoVar)==row(CoVar)]=0
  Total.SE=sqrt(sum(CoVar)+sum(VarP))
  Total.Reserve=sum(P)
  Incr=data.frame(newdata, Y, VarY, P,seP, CV=seP/P)
  out=list(Forecast=Incr,
           Totals=data.frame(Total.Reserve,
                             Total.SE=Total.SE,
                             CV=Total.SE/Total.Reserve))
  return(out)}

tail.years=6
fdat=data.frame(
  origin=rep(2007:2013, n+tail.years),
  dev=rep(1:(n+tail.years), each=n)
)
fdat=within(fdat,{
  cal=origin+dev-1
  a7=ifelse(origin==2013,1,0)
  a6=ifelse(origin==2012,1,0)
  originf=factor(origin)
  p34=ifelse(dev<2,1,0)
  d1=ifelse(dev<2,1,0)
  d27=ifelse(dev<2,0,dev-1)
})

reserve2=log.incr.predict(fit2, subset(fdat, cal>2013))
reserve2$Totals

reserve3=log.incr.predict(fit3, subset(fdat, cal>2013))
reserve3$Totals


round(xtabs(P~origin+dev, reserve3$Forecast))


round(summary(MackChainLadder(cum.triangle, est.sigma='Mack',
                              tail=1.05, tail.se=0.02))$Totals,2)


(cum.triangle.ny=t(apply(inc.triangle.ny,1,cumsum)))

f.ny=sapply((n-1):1, function(i){
  sum(cum.triangle.ny[1:(i+1), n-i+1])/sum(cum.triangle.ny[1:(i+1), n-i])
})

(f.ny=c(f.ny[-(n-1)],1))

full.triangle.ny=cum.triangle.ny
for (k in 2:(n-1)){
  full.triangle.ny[(n-k+2):n, k+1]=full.triangle.ny[(n-k+2):n,k]*f[k]
}
(sum(re.reserve.995=full.triangle.ny[,n]-rev(latest.paid)))


library(sp)
library(CASdatasets)

data(credit)
credit=credit[,-20]

table(credit$class)





# credit$class=credit$class-1

## Descripsion of the data

credit_f=credit
credit_f$age=cut(credit_f$age, c(0,25,Inf))
credit_f$credit_amount=cut(credit_f$credit_amount, c(0,4000,Inf))
credit_f$duration=cut(credit_f$duration, c(0,15,36,Inf))


# install.packages('rcompanion')

#calculate Cramer's V
# cramerV(data)
#cramer=function(i) sqrt(chisq.test
# (table(credit.f[,i], credit_f$class))$statistic/length(credit_f[,i]))


# 
# library(rcompanion)
# cramerV(credit_f)
# 크라머v correlation 

k=ncol(credit_f)-1
cramer=function(i) sqrt(chisq.test(
  table(credit_f[,i],credit_f$class))$statistic/(length(credit_f[,i])))

pv=function(i) chisq.test(table(credit_f[,i],credit_f$class))$p.value

CRAMER=data.frame(variable=names(credit)[1:k],
                  cramerv=Vectorize(cramer)(1:k),
                  p.value=Vectorize(pv)(1:k))
vCRAMER=CRAMER[order(CRAMER[,2],decreasing=TRUE),]

##상관분석 시각화
par(mar=c(10,4,4,0))
barplot(vCRAMER[,2],names.arg=vCRAMER[,1],las=3)



aggregate(credit[,c('age','duration')],by=list(class=credit$class),mean)

Q=quantile(credit$age,seq(0,1,by=.1))
Q[1]=Q[1]-1
cut.age=cut(credit$age,Q)
(prop=prop.table(table(cut.age,credit$class),1))

barplot(t(prop))
abline(h=mean(credit$class==0), lty=2)
## 분위수 시각화 Proportion of class categories, per age category



#scoring

Y=credit$class
reg=glm(class~age+duration, data=credit, family=binomial(link='probit'))
summary(reg)


S=pnorm(predict(reg))
S

FP= function(s) sum((S>s)*(Y==0))/sum(Y==0)*100
TP= function(s) sum((S>s)*(Y==1))/sum(Y==1)*100
u=seq(0,1,length=251)
plot(Vectorize(FP)(u), Vectorize(TP)(u),types='s',
     xlab='False Positive Rate (%)', ylab='True Positive Rate (%)')
abline(a=0, b=1,col='grey')


# install.packages('ROCR')
library(ROCR)

pred=prediction(S,Y)
perf=performance(pred,'tpr','fpr')
plot(perf)

# install.packages('verification')
library(verification)

roc.plot(Y,S,Xlab='False Positive Rate',
         ylab='True Positive Rate', main='', CI=TRUE,
         n.boot=100, plot='both', binormal=TRUE)

library(pROC)

roc=plot.roc(Y,S,main='', percent=TRUE, ci=TRUE)
roc.se=ci.se(roc,specificities=seq(0,100,5))
plot(roc.se, type='shape', col='grey')


#KS statistic 
max(attr(perf,'y.values')[[1]]-attr(perf,'x.values')[[1]])


plot(ecdf(S[credit$class==0]),main='',xlab='',pch=19, cex=.2)
plot(ecdf(S[credit$class==1]), pch=19,cex=.2, col='grey', add=TRUE)
legend('bottomright',c('Score=0','Score=1'), pch=19, col=c('black','grey'),lwd=1,bty='n')
perf=performance(pred,'tpr','fpr')
ks=perf@y.values[[1]]-perf@x.values[[1]]
(seuil=pred@cutoffs[[1]][which.max(ks)])
arrows(seuil,1-perf@y.values[[1]][which.max(ks)],seuil,
       1-perf@x.values[[1]][which.max(ks)],col='black',length=0.1,code=3)

performance(pred,"ecost")


performance(pred,'auc')@y.values


library(car)

credit_rcd=credit_f
credit_rcd$checking_status=
  recode(credit_rcd$checking_status,
         "'A14'='No checking account'; 'A11'='CA < 0 euros';'A12'=
         'CA in [0-200 euros[';'A13'='CA>200 euros'")
credit_rcd$credit_history=
  recode(
    credit_rcd$credit_history,
    "c('A30','A31')='critical account';
    c('A32','A33')='existing credits paid back duly till now';
    'A34'='all credits paid back duly'")
credit_rcd$purpose=recode(
  credit_rcd$purpose,"'A40'='Car (new)';
  'A41'='Car (used)';c('A42', 'A43', 'A44','A45')='Domestic equipment';c('A46','A48'
  ,'A49')='Studies-Business';
  'A47'='Holidays';else='Else'")
credit_rcd$savings=recode(
  credit_rcd$savings,"c('A65', 'A63', 'A64')=
  'No savings or >500 euros'; c('A62', 'A61')='<500 euros'")
credit_rcd$employment=recode(
  credit_rcd$employment,"c('A71','A72')='unemployed or <1 year';'A73'=
  'E[1-4[years';c('A74','A75')='>4 years'")

credit_rcd$personal_status=
  recode(credit_rcd$personal_status,
         "'A91'='male divorced/separated';'A92'='female divorced/separated/married';
         c('A93','A94')='male single/married/widowed';'A95'='female : single'")
credit_rcd$property_magnitude=
  recode(credit_rcd$property_magnitude,"'A121'='Real estate';'A124'='No property';
         else='Else'")
credit_rcd$otherpayment_plans=
  recode(credit_rcd$other_payment_plans,"'A143'='None';else='Banks-Stores'")


credit_rcd$housing=(
  recode(credit_rcd$housing,"'A152'='Owner';else='Else'")
)


# Training and Testing samples

set.seed(123)
index=sort(sample(nrow(credit), 644, replace=FALSE))
table(credit$class[index])


install.packages('rngWELL')
install.packages('randtoolbox')
library(rngWELL)
library(randtoolbox)


set.generator(name='congruRand', mod=2^(31)-1, mult=397204094, incr=0,seed=123)
U=runif(1000)
index=sort(which(rank(U)<644))
table(credit$class[index])


train_db=credit_rcd[index,]
valid_db=credit_rcd[-index,]



#Logistic Regression

# reg=glm(Y~X1+X2+X3, data=df, family=binomial(link='logit'))


Y=credit[index,'class']
X=as.matrix(cbind(1,credit[index,c('age','duration')]))

##newton-raphson
beta=as.vector(lm(Y~0+X[,1]+X[,2]+X[,3])$coefficients)
BETA=NULL
for( s in 1:6){
  pi=exp(X%*%beta)/(1+exp(X%*%beta))
  gradient=t(X)%*%(Y-pi)
  omega=matrix(0,nrow(X),nrow(X));diag(omega)=(pi*(1-pi))
  hessian=-t(X)%*%omega%*%X
  beta=beta-solve(hessian)%*%gradient
  BETA=cbind(BETA,beta)
}


BETA



(SD=sqrt(diag(solve(-hessian))))

cbind(BETA[,6],SD,BETA[,6]/SD, 2*(1-pnorm(abs(BETA[,6]/SD))))


##glm object

reg=glm(class~age+duration, data=credit[index,],
        family=binomial(link='logit'))
summary(reg)



##Logistic Regressio on Categorial variates

reg=glm(class~credit_history, data=credit, family=binomial(link='logit'))
summary(reg)



cbind(prop.table(table(credit$credit_history, credit$class),1),
      logit=predict(reg,newdata=data.frame(credit_history=
                                             levels(credit$credit_history)),
                    type='response'))



length(credit_rcd$credit_history)
length(credit_rcd$purpose)

credit_rcd$credit_history <- factor(credit_rcd$credit_history,
                                    levels=c("all credits paid back duly",
                                             "critical account",
                                             "existing credits paid back duly till now"))
credit_rcd$purpose <- factor(credit_rcd$purpose,
                             levels=c("Car (new)","Car (used)",
                                      "Domestic equipment", "Else","Studies-Business"  ))

reg=glm(class~credit_history*purpose, data=credit_rcd,
        family=binomial(link='logit'))

p.class=matrix(predict(reg,newdata=data.frame(
  credit_history=rep(levels(credit_rcd$credit_history), each=length(levels(credit_rcd$purpose))),
  purpose=rep(levels(credit_rcd$purpose), length(levels(credit_rcd$credit_history)))), type='response'),
  ncol=length(levels(credit_rcd$credit_history)),
  nrow=length(levels(credit_rcd$purpose)))

rownames(p.class)=levels(credit_rcd$purpose)
colnames(p.class)=levels(credit_rcd$credit_history)


reg=glm(class~credit_history*purpose, data=credit_rcd)

p.class.linear=matrix(predict(reg,newdata=data.frame(
  credit_history=rep(levels(credit_rcd$credit_history), each=length(levels(credit_rcd$purpose))),
  purpose=rep(levels(credit_rcd$purpose), length(levels(credit_rcd$credit_history)))), type='response'),
  ncol=length(levels(credit_rcd$credit_history)),
  nrow=length(levels(credit_rcd$purpose)))

rownames(p.class.linear)=levels(credit_rcd$purpose)
colnames(p.class.linear)=levels(credit_rcd$credit_history)


p.class.linear

p.class.linear/(1-p.class.linear)

AIC(reg, k=2)
AIC(reg,k=log(nobs(reg)))



predictors=names(credit_rcd) [-grep('class', names(credit_rcd))]
formula=as.formula(paste('y~',paste(names(credit_rcd[,predictors]), collapse='+')))


logit=glm(class~1, data=train_db, family=binomial)

selection=step(logit, direction='forward', trace=TRUE, k=log(nrow(train_db)),
               scope=list(upper=formula))



logit=glm(class~.,data=train_db[,c('class', predictors)], family=binomial)
selection=step(logit, direction='backward', trace=TRUE, k=log(nrow(train_db)))


train_db$forward.bic=predict(selection, newdata=train_db, type='response')
valid_db$forward.bic=predict(selection, newdata=valid_db, type='response')


pred_train=prediction(train_db$forward.bic,train_db$class)
pred_valid=prediction(valid_db$forward.bic, valid_db$class)
perf_train=performance(pred_train,'tpr', 'fpr')
perf_valid=performance(pred_valid, 'tpr', 'fpr')
plot(perf_train,col='grey', lty=2, main='Forward selection')
plot(perf_valid,add=TRUE)



train_db=credit_f[index,]
valid_db=credit_f[-index,]
y=as.numeric(train_db[,'class'])
x=data.frame(model.matrix(~.,data=train_db[, -which(names(train_db)=='class')]))
# install.packages("leaps")
library(leaps)
# selec=leaps(x,y,method='Cp', nbest=1, strictly.compatible=FALSE)
#계산시간 걸린다
selec$label




plot(selec$size-1, selec$Cp, xlab='Number of predictors', ylab='Cp')
selec


best.model=selec$which[which((selec$Cp==min(selec$Cp))),]
z=cbind(x,y)

formula=as.formula(paste('y~', paste(colnames(x)[best.model],
                                     collapse='+')))
logit=glm(formula, data=z, family=binomial(link='logit'))
summary(logit)

xp=data.frame(model.matrix(~.,data=valid_db[-which(names(valid_db)=='class')]))
xp
predclass=predict(logit, xp,type='response')
pred=prediction(predclass,valid_db$class,label.ordering=c(0,1))
performance(pred,'auc')@y.values[[1]]


pred.train=prediction(predict(logit, newdata=x,type='response'), train_db$class)
pred.valid=prediction(predict(logit, newdata=xp,type='response'),valid_db$class)
perf.train=performance(pred.train,'tpr','fpr')
perf.valid=performance(pred.valid,'tpr', 'fpr')
plot(perf.train, col='grey', lty=2,main='Logit model selected with leaps')
plot(perf.valid,add=TRUE)
legend('bottomright', c('Training dataset', 'Validation dataset'),
       lty=c(2,1), col=c('grey', 'black'),lwd=1)


regglm=glm(class~age+duration, data=credit[index,], family=binomial(link='logit'))
summary(regglm)

library(nlme)
library(mgcv)
# reggam=gam(Y~s(x1),data=db, family=binomial(link='logit'))
# reggam=gam(Y~s(x))
reggam=gam(class~s(age,duration),data=credit[index,],family=binomial(link='logit'))


pglm=function(x1,x2){
  return(predict(regglm,newdata=data.frame(duration=x1,age=x2),type='response'))}
pgam=function(x1,x2){
  return(predict(reggam,newdata=data.frame(duration=x1,age=x2), type='response'))}
M=31
cx1=seq(from=min(credit$duration), to=max(credit$duration),length=M)
cx2=seq(from=min(credit$age), to=max(credit$age), length=M)
Pgam=outer(cx1,cx2,pgam)
Pglm=outer(cx1,cx2,pglm)
persp(cx1,cx2,Pgam)
contour(cx2,cx2,Pgam)
persp(cx1,cx2,Pglm)
contour(cx1,cx2,Pglm)



pkmeans=function(x1,x2,k=25){
  D=as.matrix(dist(rbind(credit[index,c('duration','age')],
                         c(x1,x2))))[length(index)+1,1:length(index)]
  i=as.vector(which(D<=sort(D)[k]))
  return(mean((credit[index,'class'])[i]))}


install.packages('glmnet')
library('glmnet')

ridge=glmnet(x,y,alpha=0, family='binomial', lambda=c(0,1,2), standardize=TRUE)
train.db=credit_rcd[index,]

yt=as.numeric(train.db[,'class'])
xt=model.matrix(~.,data=train.db[,-which(names(train.db)=='class')])
set.seed(1)
cvfit=cv.glmnet(xt,yt,alpha=0, family='binomial', type='auc',nlambda=100)
cvfit$lambda.min


plot(cvfit)
abline(v=log(cvfit$lambda.min), col='blue', lty=2)

fits=glmnet(xt,yt,alpha=0,family='binomial', lambda=seq(cvfit$lambda[1],
                                                        cvfit$lambda[100], length=10000)
            ,standardize=TRUE)
plot(fits,xvar='lambda', label='T')

# train.db=credit_rcd[index,]
valid.db=credit_rcd[-index,]
yv=as.numeric(valid.db[,'class'])
xv=model.matrix(~.,data=valid.db[,-which(names(train.db)=='class')])

yvpred=predict(fits,newx=xv,type='response')
library(ROCR)
roc=function(x) {performance(prediction(yvpred[,x],yv),'auc')@y.values[[1]]}
vauc=Vectorize(roc)(1:ncol(yvpred))
fits$lambda[which.max(vauc)]

plot(fits$lambda,vauc,type='l')

vauc[which.max(vauc)]



lasso=glmnet(x,y,alpha=1,family='binomial', lambda=c(0,1,2),standardize = TRUE)

set.seed(235)
cvfit=cv.glmnet(xt,yt,alpha=1, family='binomial',type='auc', nlambda=100)
cvfit$lambda.min

plot(cvfit)
abline(v=log(cvfit$lambda.min),col='grey', lty=2)


fitx=glmnet(xt,yt,alpha=0,family='binomial', lambda=seq(cvfit$lambda[1],
                                                        cvfit$lambda[71],length=10000),
            standardize = TRUE)
plot(fits,xvar='lambda',label='T')



X1=credit[,'age']
X2=credit[,'duration']
Y=credit[,'class']
credit

criteria1=function(s){
  sum((Y[X1<=s]-mean(Y[X1<=s]))^2)+sum((Y[X1>s]-mean(Y[X1>s]))^2)}
criteria2=function(s){
  sum((Y[X2<=s]-mean(Y[X2<=s]))^2)+sum((Y[X2>s]-mean(Y[X2>s]))^2)}

S=seq(0,100,length=501)
plot(S,Vectorize(criteria2)(S), type='l', ylab='Sum of squares',xlab='Splitting point')
lines(S,Vectorize(criteria1)(S),type='l',lty=2)
legend(70,205,c('Variable X1', 'Variable X2'),lty=2:1)


S[which(Vectorize(criteria2)(S)==min(Vectorize(criteria2)(S)))]


s.star=34.5
creteria1.lower=function(s){
  sum((Y[(X1<=s)&(X2<=s.star)]-mean[((X1<=s)&(X2<=s.star)]))^2)}

  
  
  


# install.packages('forecast')
# install.packages('demography')
# install.packages('gnm')
# install.packages('StMoMo')
# install.packages('fanplot')
library(forecast)
library(demography)
library(gnm)
library(StMoMo)
library(fanplot)
library(ggplot2)


constLC <- function(ax,bx,kt, b0x, gc, wxt, ages){
  c1 <- mean(kt[1,],na.rm=TRUE)
  c2 <- sum(bx[,1], na.rm=TRUE)
  list(ax <- ax+c1*bx[,1],
       bx[,1] <- c2*(kt[1,]-c1))
}


USdata <- hmd.mx(country='USA', username='donmany85@gmail.com', password='******!')



LC <- StMoMo(link='log', staticAgeFun=TRUE,
             periodAgeFun='NP', constFun=constLC)
LC <- lc()

Ext <- USdata$pop$male 
Dxt <- USdata$rate$male*Ext 
ages <- USdata$age 
years <- USdata$year

ages.fit <- 0:110
years.fit <- 1933:2010
LCfit <- fit(LC, Dxt=Dxt, Ext=Ext, ages=ages, years=years,
             ages.fit=ages.fit, years.fit=years.fit)
plot(LCfit)
prod1 <- as.data.frame(LCfit$bx%*%LCfit$kt)
LCfitax <- as.data.frame(LCfit$ax)
LCfitaxmatrix <- as.data.frame(replicate(78,LCfitax))
LCfitdeathrates <- LCfitaxmatrix+prod1




# USdrates1940 <- as.data.frame(USdeathratestot[1:111,8])

LCfitdeathrates1940 <- as.data.frame(LCfitdeathrates[1:111,8])

ages1940 <- 0:110 
dim(ages1940)

plot(LCfitdeathrates[1:111,10], type='o', lty=3, col='red')


CBD <- cbd(link='log')
CBDfit <- fit(CBD, Dxt=Dxt, Ext=Ext, ages=ages, years=years,
              ages.fit=ages.fit, years.fit=years.fit)


plot(CBDfit)

APC <- apc()
APCfit <- fit(APC, Dxt=Dxt, Ext=Ext, ages=ages, years=years,
              ages.fit=ages.fit, years.fit=years.fit)



M7 <- m7(link='log')
M7fit <- fit(M7, Dxt=Dxt, Ext=Ext, ages=ages, years=years,
             ages.fit=ages.fit, years.fit=years.fit)


AIC(LCfit)
BIC(LCfit)


AIC(CBDfit)
BIC(CBDfit)


LCres <- residuals(LCfit)

plot(LCres, type='colourmap', main='LC residuals')
plot(LCres)


LCfor <- forecast(LCfit, h=50)
plot(LCfor$fitted)

LCsim <- simulate(LCfit, nsim=1000, h=50)
plot(LCfit$years, LCfit$kt[1,], type='l', xlim=c(1933,2060), ylim=c(-100,50))
matlines(LCsim$kt.s$years, LCsim$kt.s$sim[1,1:50,1:50], type='l')


mxt <- LCfit$Dxt/LCfit$Ext
(plot1 <- plot(LCfit$years, mxt['65',], type='l'))
(plot2 <- plot(LCfit$years, mxt['75',], type='l'))
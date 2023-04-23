

# install.packages("CASdatasets",repos="http://dutangc.free.fr/pub/RRepos/", type="source")
library(CASdatasets)
data(freMTPL2freq)

attach(freMTPL2freq)

summary(freMTPL2freq)
summary(ClaimNb)
hist(ClaimNb)

hist(Density)

freMTPL2freq['VehAgeGrp'] <- cut(freMTPL2freq$VehAge, c(0,1,11,Inf),
                                 include.lowest=TRUE, right=FALSE)

freMTPL2freq['DrivAgeGrp'] <- cut(freMTPL2freq$DrivAge,
                                  c(18,21,26,31,41,51,71,100),
                                  include.lowest=TRUE, right=FALSE)

set.seed(100)

ll <- sample(c(1:nrow(freMTPL2freq)), round(0.9*nrow(freMTPL2freq)), replace=FALSE)
learn <- freMTPL2freq[ll,]
test <- freMTPL2freq[-ll,]

ll

frequ1 <- formula(learn$ClaimNb~learn$VehPower+learn$VehAgeGrp+
                    learn$DrivAgeGrp+learn$BonusMalus+learn$Density+learn$VehBrand+
                    learn$VehGas+learn$Area+learn$Region+offset(log(learn$Exposure)))
glm1 <- glm(frequ1, data=learn, family=poisson())
summary(glm1)


frequ2 <- formula(learn$ClaimNb~learn$VehPower+learn$VehAgeGrp+
                    learn$DrivAgeGrp+learn$BonusMalus+learn$Density+
                    learn$VehGas+learn$Region+offset(log(learn$Exposure)))
glm2 <- glm(frequ2, data=learn, family=poisson())

summary(glm2)


frequ3 <- formula(learn$ClaimNb~learn$VehPower+learn$VehAgeGrp+
                   learn$DrivAgeGrp+learn$BonusMalus+learn$Density+
                   learn$VehGas+learn$Area+learn$Region+offset(log(learn$Exposure)))
glm3 <- glm(frequ3, data=learn, family=poisson())
summary(glm3)


# CValidation
learn$fit <- fitted(glm1)
test$fit <- predict(glm1, newdata=test, type='response')

inSampleLoss <- 2*(sum(learn$fit)-sum(learn$ClaimNb)
                   +  sum(log((learn$ClaimNb/learn$fit)^(learn$ClaimNb))))
inSampleLoss


OutOfSampleLoss <- 2*(sum(test$fit)-sum(test$ClaimNb)
                      +sum(log((test$ClaimNb/test$fit)^(test$ClaimNb))))
OutOfSampleLoss


library(rpart)
# install.packages("rpart.plot")
library(rpart.plot)

tree <- rpart(cbind(Exposure, ClaimNb)~VehPower+
                VehAgeGrp+DrivAgeGrp+BonusMalus+
                Density+VehBrand+VehGas,
              data=learn, method='poisson',
              control=rpart.control(rval=1,
                                    minbucket=7000,
                                    cp=0.0005))

rpart.plot(tree) 
summary(tree)




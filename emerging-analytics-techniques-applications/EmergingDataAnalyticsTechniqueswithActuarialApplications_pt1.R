library(tidyverse)
library(ChainLadder)

cs1 <- read.csv(file='./MedicalCare.csv', header=TRUE)

cs1_tbl <- aggregate(Payments~Month+Delay, data=cs1,sum) %>% spread(key='Delay', value='Payments')
cs1_tbl2 <- cs1_tbl[24:nrow(cs1_tbl),]
rownames(cs1_tbl2) <- cs1_tbl2[,1]
cs1_tbl2 <- cs1_tbl2[,-1]

cs1_tri <- as.triangle(as.matrix(cs1_tbl2), origin='origin', dev='dev', value='value')
cs1_tri_cum <- incr2cum(cs1_tri)

plot(cs1_tri_cum, lattice=FALSE, ylab='Cumulative Payments')
plot(cs1_tri_cum, lattice=TRUE, ylab='Cumulative Payments(TRUE-lattice)')

g <- attr(ata(cs1_tri_cum), 'vwtd')
g <- c(g,1)
full_cs1 <- cbind(cs1_tri_cum, Ult=rep(0,13))

n <- nrow(full_cs1)
for(k in 1:n){
  full_cs1[(n-k+1):n, k+1] <- full_cs1[(n-k+1):n,k]*g[k]
}


sum(full_cs1[,14]-getLatestCumulative(cs1_tri_cum))

Pd_to_Dt <- getLatestCumulative(cs1_tri_cum)

linkratios <- c(attr(ata(cs1_tri_cum), 'vwtd'), tail=1.000)
round(linkratios,3)
LDF <- rev(cumprod(rev(linkratios)))
names(LDF) <- colnames(cs1_tbl2)
round(LDF,3)


EstUlt <- Pd_to_Dt*rev(LDF)
Reserve <- EstUlt-Pd_to_Dt
exhibit <- data.frame(Pd_to_Dt,LDF=round(rev(LDF),3),EstUlt, Reserve)
exhibit <- rbind(exhibit, data.frame(Pd_to_Dt=sum(Pd_to_Dt), LDF=NA, EstUlt=sum(EstUlt), Reserve=sum(Reserve),
                                     row.names='Total'))
exhibit 


mack1 <- MackChainLadder(cs1_tri_cum,est.sigma='Mack')
plot(mack1, lattice=FALSE)
mack1


# FFT (f, M, param, distN)
# input: f(sezuence), M(truncation point), distN( the dist. of N)
# and param(list with the parameters of the severity dist.)
# output: vector with the the corresponding single probabilities
# q[k]=P(S_N=k)


FFt <- function(f,M, param, distN){
  fhat <- fft(f, inverse=FALSE)
  
  #apply the generating function of N to fhat
  if(distN=='Poi'){
    #param <- c(lambda)
    qhat <- exp(param[1]*(fhat-1))
  }
  if(distN=='Bin'){
    #param <- c(n,p)
    qhat <- (1-param[2]+param[2]*fhat)^param[1]
  }
  if(distN=='NB'){
    #param <- c(r,p)
    qhat <- (param[2]/(1-(1-param[2])*fhat))^param[1]
  }
  
  q <- 1/M*fft(qhat, inverse=TRUE)
  return(q)
}



# FFT.tilting(f, M, param, theta, distN)
# 
# input: f(sequence), M(trancationn point), distN(the dist. of N),
# param(list with the parameters of the severity dist.)
# 
# output: vector with the the corresponding single probabilities
# q[k]=P(s_N=k)


FFT.tilting <- function(f, M, param, theta, distN){
  fnew <- exp(-theta*(0:(M-1)))*f
  fhat <- fft(fnew, inverse=FALSE)
  
  
  #apply the generating function of N to fhat  
  if(distN=='Poi'){
    #param <- c(lambda)
    qhat <- exp(param[1]*(fhat-1))
  }
  if(distN=='Bin'){
    #param <- c(n,p)
    qhat <- (1-param[2]+param[2]*fhat)^param[1]
  }
  if(distN=='NB'){
    #param <- c(r,p)
    qhat <- (param[2]/(1-(1-param[2])*fhat))^param[1]
  }
  q <- exp(theta*(0:(M-1)))*(1/M*fft(qhat, inverse=TRUE))
  return(q)}











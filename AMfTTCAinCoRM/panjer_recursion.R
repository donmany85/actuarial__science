# input:lambda(parameter for the Poisson distribution)
# output:vector with the corresponding single probabilities
# q[k]=P(S_N=k)


Panjer.Poisson <- function(lambda, f){
  #check, if necessary requirements are satisfied
  if(sum(f) !=1 |any(f<0)){
    stop('f is not a density.')
  }
  if (lamba*sum(f)>727){
    stop('underflow.')
  }
  if(length (f)==3&f[2]==0&f[3]==1){
    stop('P(SN=0)=1')
  }
  
  s <- rep(0,100000)
  k <- 1
  
  #first term of the recurssion
  cumul <- q <- f0 <- exp(-lambda*(1-f[1]))
  l <- length(f)
  #run the calculation until q[0]+q[1]+....q[k]>0.99999 holds
  while(cumul<0.99999){
    for (j in 1:min(l-1,k)){
      s[j] <- lambda/k*j*f[j+1]*q[k+1-j]
    }
    q <- c(q,sum(s))
    cumul <- cumul+sum(s)
    k <- k+1
    }
  return(q)}


# input: n, p(parameter for the Binomial Dist.) and f(severity dist.)
# vector with the corresponding single probabilities 
# q[k]=P(S_N=k)





Panjer.Bin <- function(n,p,f){
  #check, if neccesary requirements are satisfied
  if (sum(f) !=1|any(f<0)){
    stop('f is not a density.')
  }
  if (length(f)==3&f[2]==0&f[3]==1){
    stop('P(SN=0)=1')
  }
  if(p>=1|p<=0){
    stop('p is not a probability.')
  }
  
  s <- rep(0,100000)
  k <- 1
  a <- p/(p-1)
  b <- p*(n+1)/(1-p)
  
  #first term of the recursion
  cumul <- q <- f0 <- ((a-1)/(a*f[1]-1))^((a+b)/a)
  l <- length(f)
  
  # run the calculation until q[0]+q[1]+...+q[k]>0.99999 holds
  while(cumul<0.9){
    for (j in 1:min(l-1,k)){
      s[j] <- (1/(1-a*f0))*(a+j*b/k)*f[j+1]*q[k+1-j]
    }
    q <- c(q,sum(s))
    cumul <- cumul+sum(s)
    k <- k+1
  }
  return(q)}


# dirichlet nils example

library(fastperm)
n = 1000

set.seed(102)

ts = seq(-0.102,-1.02/2, length.out = n)

U = runif(n)

I =  U >= 1/3 * pnorm((ts+2)/0.7) + 2/3 * pnorm((ts-1)/0.7)

xleft = ts[I]
xright = ts[!I]

nleft = length(xleft)
nright = length(xright)

a = c(rep(-Inf, nright), xleft)
b = c(xright, rep(Inf, nleft))


T = 100

alpha_param = 1

X = matrix(NA, nrow = T, ncol = n)
# generate samples from prior
for (t in 1:T) {
  X[t,1] = rnorm(1)
  Us = runif(n-1)
  for (i in 2:n) {
    if(Us[i-1] <= alpha_param / (alpha_param +i-1)){
      X[t,i] = rnorm(1)
      
      
    }else{
      if(i==2){
        X[t,i] = X[t,i-1]
      }else{
        X[t,i] = sample(x = X[t,1:(i-1)],1)
      }
      if(X[t,i]==1){
        print("her")
      }
    }
    
  }
}


nonzeroperms = rep(TRUE,T)
for (t in 1:T) {
  for (i in 1:n){
    if(X[t,i]< a[i] || X[t,i] > b[i]){
      nonzeroperms[t] = FALSE
      #print(i)
    }
  }
  
}


# nils has now generated T marginal samples from the *marginal* Dirichlet preusess

logperms = get_log_permanents(X,a,b, FALSE)



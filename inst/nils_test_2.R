# this is my *nils test* osv
library(Fastperm)
filename = "/Users/peraugust/Downloads/master_debug.csv"

n = 200
dat = read.csv(filename)
dat

as = dat[,2:(n+1)]
bs = dat[,(n+2):(2*n+1)]
xs = dat[,(2*n+2):(3*n+1)]
alphas = dat[,(3*n+2):(4*n+1)]
betas = dat[,(4*n+2):(5*n+1)]
gammas = dat[,(5*n+2):(6*n+1)]

T = 988
diffs = rep(0, T)

for (t in 1:T) {
  #print(t)
  nils = test_get_alphabetagamma_input(xs[t,],as[t,],bs[t,],debug=FALSE)
  # print(sum(abs(nils$alpha - alphas[t,])))
  # print(sum(abs(nils$beta - betas[t,])))
  # print(sum(abs(nils$gamma- gammas[t,])))
  
  ss = sum(abs(nils$alpha - alphas[t,])) + sum(abs(nils$beta - betas[t,]))+ sum(abs(nils$gamma- gammas[t,]))
  k_true = sum(alphas[t,]>0)
  if(ss>0){
    print(t)
  }
  if(k_true != nils$k){
    print("Error at k, true and given:")
    print(k_true)
    print(nils$k)
  }
}
diffs
max(diffs)
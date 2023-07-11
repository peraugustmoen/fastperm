# this is my *nils test* osv
library(Fastperm)
filename = "/Users/peraugust/Downloads/halla.csv"

n = 8
dat = read.csv(filename)
dat

as = dat[,2:(n+1)]
bs = dat[,(n+2):(2*n+1)]
xs = dat[,(2*n+2):(3*n+1)]
logperms = dat[,3*n+2]

T = 1000
diffs = rep(0, T)

for (t in 1:T) {
  print(t)
  nils = get_log_permanents(as.numeric(xs[t,]),as.numeric(as[t,]),as.numeric(bs[t,]), FALSE)
  print(abs(nils - logperms[t]))
  diffs[t]  = abs(nils - logperms[t])
}
diffs
max(diffs)

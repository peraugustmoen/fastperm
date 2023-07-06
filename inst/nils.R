
library(Fastperm)
n = 5
Tt = 200
points = seq(-3,3, length.out = n+5)

a = points[1:n]
b = points[6:(n+5)]

a[1:round(n/2)] = -Inf
b[(round(n/2)+1):n] = Inf

x = rnorm(n*Tt)

X = matrix(x, ncol = n)
#X = matrix(c(1.9611550 , 0.3144615 , 0.8194104 , 1.0274892, -1.7606855),ncol=n)

print("a:")
print(a)
print("b:")
print(b)

a = get_log_permanents(X,a,b, TRUE)
a



a= c(-Inf ,-Inf ,-1.666667 ,-1.000000 ,-0.333333 )
b =c(0.333333 ,1.000000 ,Inf ,Inf ,Inf )
x = c(-0.952536, -0.496228 ,0.825848 ,0.913911 ,0.925020 )
aa = get_log_permanents(X,a,b, TRUE)

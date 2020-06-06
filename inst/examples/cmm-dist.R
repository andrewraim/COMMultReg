library(COMMultReg)

p = c(1,2,3) / 6
y = r_cmm(n = 100, m = 10, p = p, nu = 0.8, burn = 1000, thin = 10)
head(y)

library(COMMultReg)
library(MM)

# ----- Read the pollen dataset -----
data("pollen")
pollen = as.data.frame(pollen)
colnames(pollen) = c("Pine", "Fir", "Oak", "Alder")
head(pollen)

y = as.matrix(pollen)
n = nrow(pollen)
k = ncol(pollen)
X = matrix(1, n, 1)
W = matrix(1, n, 1)

# ----- Fit CMM distribution -----
ctrl = cmm_reg_control(verbose = TRUE, base = 1)
out_cmm = cmm_reg(cbind(Pine,Fir,Oak,Alder) ~ 1, data = pollen, control = ctrl)

# Alternatively, use the "raw" interface
# out_raw_cmm = cmm_reg_raw(y = y, X = X, W = W, control = ctrl)

# Print summary
print(out_cmm)

# Get fitted values
fit_out = fitted(out_cmm, newX = X, newW = W)
p_hat = as.numeric(fit_out[1, 1:k])
nu_hat = fit_out[1, k+1]

# Estimates of coefficient and their covariances
psi_hat = coef(out_cmm)
V_hat = vcov(out_cmm)

# Log-likelihood
logLik(out_cmm)

# E(Y) and Var(Y) under the estimated model, for one observation
e_cmm(as.integer(m[1]), p_hat, nu_hat)
v_cmm(as.integer(m[1]), p_hat, nu_hat)

# Compare to empirical mean and variance
colMeans(y)
var(y)

library(COMMultReg)
library(reshape2)
library(vcdExtra)

data(Alligator)
gator = dcast(Alligator, lake + sex + size ~ food, sum, value.var = "count")

y = as.matrix(gator[,c("fish", "invert", "reptile", "bird", "other")])
m = rowSums(y)
k = ncol(y)
n = nrow(y)

# ----- Try some model selection on W with CMM -----
X = model.matrix(~ size + lake, data = gator)
W_levels = list(
	model.matrix(~ 1, data = gator),
	# One covariate
	model.matrix(~ sex, data = gator),
	model.matrix(~ size, data = gator),
	model.matrix(~ lake, data = gator),
	# Two covariates
	model.matrix(~ sex + size, data = gator),
	model.matrix(~ sex + lake, data = gator),
	model.matrix(~ size + lake, data = gator),
	# Three covariates
	model.matrix(~ sex + size + lake, data = gator)
)
aic_levels = numeric(length(W_levels))

ctrl = cmm_poisreg_control(verbose = FALSE)
for (l in 1:length(W_levels)) {
	logger("Fitting model %d of %d\n", l, length(W_levels))
	W = W_levels[[l]]
	out_cmm = cmm_poisreg(y, m, X, W, base = 1, control = ctrl)
	aic_levels[l] = AIC(out_cmm)
}

# ----- Look at best fitting model -----
idx = which.min(aic_levels)
W = W_levels[[idx]]
X_names = c(
	sprintf("invert:%s", colnames(X)),
	sprintf("reptile:%s", colnames(X)),
	sprintf("bird:%s", colnames(X)),
	sprintf("other:%s", colnames(X))
)
W_names = sprintf("W:%s", colnames(W))

ctrl = cmm_poisreg_control(verbose = FALSE)
out_cmm = cmm_poisreg(y, m, X, W, base = 1,	control = ctrl)
print(out_cmm)

# Get fitted values
fit_out = fitted(out_cmm, newX = X, newW = W)
p_hat = fit_out[,-(k+1)]
nu_hat = fit_out[,k+1]

ecmm_out = matrix(NA, n, k)
colnames(ecmm_out) = sprintf("E(%s)", colnames(y))
for (i in 1:n) {
	ecmm_out[i,] = ecmm(as.integer(m[i]), as.numeric(p_hat[i,]), as.numeric(nu_hat[i]))
}

# Print observed vs expected counts under estimated CMM model
cbind(y, round(ecmm_out,4))

# GOF statistic
sum((y - ecmm_out)^2 / ecmm_out)

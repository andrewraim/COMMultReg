library(COMMultReg)
library(reshape2)
library(vcdExtra)

data(Alligator)
gator = dcast(Alligator, lake + sex + size ~ food, sum, value.var = "count")

y = as.matrix(gator[,c("fish", "invert", "reptile", "bird", "other")])
k = ncol(y)
n = nrow(y)
m = rowSums(y)

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

ctrl = cmm_reg_control(verbose = FALSE, base = 1)
for (l in 1:length(W_levels)) {
	logger("Fitting model %d of %d\n", l, length(W_levels))
	W = W_levels[[l]]
	out_cmm = cmm_reg_raw(y, X, W, control = ctrl)
	aic_levels[l] = AIC(out_cmm)
}

# ----- Look at best fitting model -----
idx = which.min(aic_levels)
W = W_levels[[idx]]

ctrl = cmm_reg_control(verbose = FALSE, base = 1)
out_cmm = cmm_reg_raw(y, X, W, control = ctrl)
print(out_cmm)

# Get fitted values
fit_out = fitted(out_cmm, newX = X, newW = W)
p_hat = fit_out[,-(k+1)]
nu_hat = fit_out[,k+1]

ecmm_out = matrix(NA, n, k)
colnames(ecmm_out) = sprintf("E(%s)", colnames(y))
for (i in 1:n) {
	ecmm_out[i,] = e_cmm(as.integer(m[i]), as.numeric(p_hat[i,]), as.numeric(nu_hat[i]))
}

# Print observed vs expected counts under estimated CMM model
cbind(y, round(ecmm_out,4))

# GOF statistic
sum((y - ecmm_out)^2 / ecmm_out)

# ----- Try the formula interface -----
ctrl = cmm_reg_control(verbose = TRUE, base = 1)
out_cmm = cmm_reg(
	formula_x = cbind(fish,invert,reptile,bird,other) ~ size + lake + sex, 
	formula_w = ~ size, 
	data = gator, control = ctrl)
print(out_cmm)

par_hat = out_cmm$par
data_xform = transform_data(out_cmm$y, out_cmm$X, out_cmm$W)
extended_intercepts(data_xform, base = 1, par_hat$beta, par_hat$gamma)
par_hat$mu

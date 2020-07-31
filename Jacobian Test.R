# write toy TMB model to check that have Jacobian right

library(TMB)
library(pscl) # use this package for nice inverse gamma functions
library(truncnorm)

dyn.unload(dynlib("SigmaTest"))
TMB::compile("SigmaTest.cpp")
dyn.load(dynlib("SigmaTest"))

Gam_Dist <- 1

data <- list(Gam_Dist=Gam_Dist)
param <- list(norm_sigma_bounds=1, norm_logsigma_Jacobian=0, 
              IG_tau_bounds=1, IG_logsigma_Jacobian=1)

obj <- MakeADFun(data, param, DLL="SigmaTest" )
obj$fn()
summary(sdreport(obj))

library(tmbstan)
library(shinystan)

fitmcmc <- tmbstan(obj, chains=3, iter=200000, thin=10,
                   lower=c(0,-Inf,0,-Inf), upper=c(Inf, Inf, Inf,Inf))
# pull out posterior vals
post <- as.data.frame(fitmcmc)

# Check truncated Normal distributions
par(mfrow=c(2,1))
x <- seq(1e-5,10, len=500)
bb <- 50
hist(post$norm_sigma_bounds, freq=FALSE, breaks=bb, 
     main = "Trunc. Normal; External Bounds",
     xlab = "SD Prior")
lines(x, dtruncnorm(x, 0,Inf))
hist(exp(post$norm_logsigma_Jacobian), freq=FALSE, breaks=bb, 
     main = "Trunc. Normal; Jacobian Adj.",
     xlab = "SD Prior")
lines(x, dtruncnorm(x, 0,Inf))

# Now Check inverse gamma priors
par(mfrow=c(2,2))
# First check that tau is gamma distributed
hist(post$IG_tau_bounds, freq=FALSE, breaks=bb, xlim = c(0, 10),
     main = "Gamma on Precision; Bounds", xlab = "Precision")
lines(x, dgamma(x,Gam_Dist,Gam_Dist))
hist(1/exp(post$IG_logsigma_Jacobia)^2, freq=FALSE, breaks=bb, xlim = c(0, 10),
     main = "Gamma on Precision; Jacobian Adj.", xlab = "Precision")
lines(x, dgamma(x, Gam_Dist,Gam_Dist))

# Also plot hist of variance, which should be inverse gamma dist., as double check
# This has some very large values that make looking at the histogram difficult
# remove very large values to help
VarPost_Bounds <- 1/post$IG_tau_bounds
VarPost_Bounds_plot <- VarPost_Bounds[VarPost_Bounds<100]
VarPost_Jacobian <- exp(2*post$IG_logsigma_Jacobian)
VarPost_Jacobian_plot <- VarPost_Jacobian[VarPost_Jacobian<100]
# need more breaks for this plot
bb <- 500
# cut off hists at 10 so can see distribution more clearly
hist(VarPost_Bounds_plot, freq=FALSE, breaks=bb, xlim = c(0, 10),
     main = "InvGamma on Variance; Bounds", xlab = "Variance")
lines(x, densigamma(x, Gam_Dist,Gam_Dist))
hist(VarPost_Jacobian_plot, freq=FALSE, breaks=bb, xlim = c(0, 10),
     main = "InvGamma on Variance; Jacobian Adj.", xlab = "Variance")
lines(x, densigamma(x, Gam_Dist,Gam_Dist))

# simulate gamma on tau to look at dist and median
tau_sim <- rgamma(100000, shape = Gam_Dist, rate = Gam_Dist)
var <- 1/tau_sim
median(var)
median(VarPost_Bounds)
median(VarPost_Jacobian)
# looks good!

# simulate inverse gamma directly on variance
var_sim <- rigamma(100000, alpha = Gam_Dist, beta = Gam_Dist)
median(var_sim) # same as above
median(VarPost_Bounds)
median(VarPost_Jacobian)
# looks good!







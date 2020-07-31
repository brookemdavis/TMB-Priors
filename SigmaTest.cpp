#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Testing two priors on variance:
  // 1) Truncated Normal(0,1) prior on the standard deviation.
  // 2) Inverse gamma prior on variance
  
  // Implement these with and without using external bounds. Without requires
  // exponentiating here and then adding a Jacobian
  // adjustment. That is done automatically when putting bounds
  // in R.
  
  DATA_SCALAR(Gam_Dist); // parameters of inverse gamma distribution
  
   // use external bounds (no need for fitting in log-space)
  PARAMETER(norm_sigma_bounds);    
  // no external bounds, need to fit in log-space to keep positive
  PARAMETER(norm_logsigma_Jacobian); 
  
  // Same with Inverse Gamma prior
  PARAMETER(IG_tau_bounds);	// use external bounds
  PARAMETER(IG_logsigma_Jacobian);	// no external bounds
  
  // set up negative log-likelihood
  Type ans = 0;

  //Truncated normal prior on standard deviation
  
  // Using external bounds (estimating sigma)
  ans -= dnorm(norm_sigma_bounds, Type(0.0), Type(1.0), true);
  // Now version without bounds, using Jacobian (estimating logsigma)
  Type norm_sigma_Jacobian = exp(norm_logsigma_Jacobian);
  ans -= dnorm(norm_sigma_Jacobian, Type(0.0), Type(1.0), true);
  // Jacobian adjustment for norm_logsigma_Jacobian
  ans -= norm_logsigma_Jacobian;
  
  // gamma distribution on tau == inverse gamma on variance
  // note that gamma in TMB is shape & scale (not shape & rate, like in JAGS)
  ans -= dgamma(IG_tau_bounds, Gam_Dist, 1/Gam_Dist, true);
  // Tau = 1/(exp(logSigma))^2
  Type IG_Tau_Jacobian = pow(exp(IG_logsigma_Jacobian), -2);
  // Now version with Jocobian
  ans -= dgamma(IG_Tau_Jacobian, Gam_Dist, 1/Gam_Dist, true);
  // Jacobian adjustment
  ans -= log(2) - 2*IG_logsigma_Jacobian;
 
  return ans;
}

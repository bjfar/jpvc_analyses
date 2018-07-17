"""Z invisible width likelihood

   Measured at LEP: we treat this as a Gaussian likelihood,
   with theory errors treated as a constrained systematic
   parameter.
"""

      const triplet<double> gamma_nu = *Dep::Z_gamma_nu;
      const triplet<double> gamma_chi_0 = *Dep::Z_gamma_chi_0;
      const double gamma_inv = gamma_chi_0.central + gamma_nu.central;
      // Average + and - errors
      const double tau_nu = 0.5 * (gamma_nu.upper + gamma_nu.lower);
      const double tau_chi_0 = 0.5 * (gamma_chi_0.upper + gamma_chi_0.lower);
      // Add theory errors in quadrature
      const double tau = std::sqrt(pow(tau_nu, 2) + pow(tau_chi_0, 2));
      lnL = Stats::gaussian_loglikelihood(gamma_inv, SM_Z::gamma_inv.mu,
        tau, SM_Z::gamma_inv.sigma, false);

name = "Z_invisible_width" 

# TODO: Hmm need to think about this, I guess the mu=0 case is for gamma_chi_0 = 0,
# i.e. mu should only scale gamma_chi_0. gamma_nu should be the background-only
# value, although I guess it varies with SM nuisance parameters. So probably
# need to set gamma_nu via "signal" parameters, but make sure it doesn't scale
# with mu, so that mu=0 only turns off the BSM contribution.

# LEP measurement data from DecayBit/include/gambit/DecayBit/SM_Z.hpp (from PDG)
gamma_inv_mu = 499.0e-3;
gamma_inv_sigma = 1.5e-3;

# Ok now build the probabilistic model for the MLE
def pars(gamma_inv,err):
    return {"loc": gamma_inv+err, "scale": gamma_inv_sigma}

def pars_nuis(err, sigma_err):
    return {"loc": err, "scale": sigma_err}

joint = jtd.JointDist([jtd.TransDist(sps.norm,pars), jtd.TransDist(sps.norm,pars_nuis)])

def get_seeds_full(samples,signal):
   Inv = samples[...,0] # invisible width measurements
   X   = samples[...,1] # Nuisance measurements (well, theory pseudo-measurements)
   return {'gamma_inv': Inv-X, 'err': X}

def get_seeds_null(samples,signal):
   Inv = samples[...,0]
   X   = samples[...,1]
   sigma_err = signal["sigma_err"]
   gamma_inv = signal["gamma_inv"]
   quad = sigma_err**2 + gamma_inv_sigma**2
   err_MLE = (1/quad) * (X*gamma_inv_sigma**2 + (Inv - gamma_inv)*sigma_err**2)
   return {'err': err_MLE}

def get_asimov(mu,signal=None):
   # Need to return data for which mu=1 or mu=0 is the MLE
   if signal is None:
      nA = 0
   else:
      nA = mu*signal['gamma_inv']
   nuis_MLEs = {'err': 0} # No nuisance parameters
   return nA, nuis_MLEs
 
nuis_options = {} # None, no nuisance fit necessary
general_options = {'BF': 0, 'error_BF': sigma} # No real need for this either since seeds give exact MLEs already.



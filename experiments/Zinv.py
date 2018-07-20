"""Z invisible width likelihood

   Measured at LEP: we treat this as a Gaussian likelihood,
   with theory errors treated as a constrained systematic
   parameter.
"""

import numpy as np
import scipy.stats as sps
import JMCtools.distributions as jtd
from JMCtools.experiment import Experiment

name = "Z_invisible_width" 

# Note: mu=0 case is for e.g. gamma_chi_0 = 0, (gamma_inv_BSM=0)
# i.e. mu should only scale gamma_chi_0. gamma_nu should be the background-only
# value, although I guess it varies with SM nuisance parameters. So
# need to set gamma_nu via "signal" parameters, but make sure it doesn't scale
# with mu, so that mu=0 only turns off the BSM contribution.

# LEP measurement data from DecayBit/include/gambit/DecayBit/SM_Z.hpp (from PDG)
gamma_inv_mu = 499.0e-3;
gamma_inv_sigma = 1.5e-3;

# We can avoid bothering with even the seed-based analytic profiling of nuisances by
# just computing the form of the profile likelihood and using that instead. Turns
# out it is proportional to a gaussian

def parfunc(gamma_inv_BSM,gamma_inv_SM,sigma_err):
    # For pseudodata generation we will always set the nuisance parameter to zero
    return {"mean": [gamma_inv_BSM+gamma_inv_SM, 0], "cov": [gamma_inv_sigma**2, sigma_err**2]}

def prof_loglike(Z,mean,cov):
    X = Z[...,0] - Z[...,1] # Should be two components, second is the nuisance parameter measurement
    return sps.norm.logpdf(X,loc=mean[0],scale=np.sqrt(np.sum(cov))) # Proportional only! Normalisation is wrong but will cancel in likelihood ratios

# Build distribution function object
mynorm = jtd.TransDist(sps.multivariate_normal) # Null transformation, just to build object
mynorm.set_logpdf(prof_loglike) # replace pdf calculation with profiled version
# Now build the joint pdf object
joint = jtd.JointDist([(jtd.TransDist(mynorm,parfunc),2)]) # make sure to let JointDist know that this is a multivariate distribution (2)

def get_seeds_full(samples,signal):
   print("samples.shape:", samples.shape)
   Inv = samples[...,0] # invisible width measurements
   X   = samples[...,1] # Nuisance measurements (well, theory pseudo-measurements)
   gamma_inv_SM = signal["gamma_inv_SM"]
   return {'gamma_inv_BSM': Inv-X-gamma_inv_SM} # Assumes that gamma_inv_SM is fixed! Which it should be.

def get_seeds_null(samples,signal):
   return {} # Nuisance parameters analytically profiled, so they don't exist as far as minimisation knows

def get_asimov(mu,signal=None):
   # Need to return data for which mu=1 or mu=0 is the MLE
   if signal is None:
      raise ValueError("'signal' dict cannot be None for Z invisible width likelihood! Must at least specific the SM contribution ('gamma_inv_SM')")
   try:
      gamma_inv_BSM = signal['gamma_inv_BSM']
   except KeyError:
      gamma_inv_BSM = 0      
   gamma_inv_SM = signal['gamma_inv_SM']
   nA = [mu*(gamma_inv_BSM) + gamma_inv_SM, 0]
   nuis_MLEs = {} # No nuisance parameters
   return nA, nuis_MLEs
 
nuis_options = {'fix_gamma_inv_SM': True, 'fix_sigma_err': True } # Make sure SM stuff stays fixed in e.g. gof test
general_options = {'gamma_inv_BSM': 0, 'error_gamma_inv_BSM': gamma_inv_sigma, **nuis_options}

# Setup done, now define the pdf and Experiment object

observed_data = np.array([gamma_inv_mu,0]) # don't forget nuisance observation! Nominally zero, by defintion. 

# Define the experiment object and options for fitting during statistical tests
e = Experiment(name,joint,observed_data,DOF=1)
 
e.define_gof_test(null_options=nuis_options,
                  full_options=general_options,
                  null_seeds=(get_seeds_null, True), # extra flag indicates that seeds are exact
                  full_seeds=(get_seeds_full, True),
                  diagnostics=None
                  )

e.define_mu_test(null_options=nuis_options,
                 null_seeds=(get_seeds_null,True),
                 scale_with_mu=['gamma_inv_BSM'],
                 )

e.define_musb_test(null_options=nuis_options,
                   mu1_seeds=(get_seeds_null,True),
                   mu0_seeds=(get_seeds_null,True),
                   scale_with_mu=['gamma_inv_BSM'],
                   asimov=get_asimov
                   )

experiments = [e]

"""Simple gaussian constraints that appear in MSSMEW.yaml"""

#gaussian_loglikelihood: obs=173.34, mu=175.04285, sigma=0.76, sigma_th=0
#  Likelihood contribution from PrecisionBit::lnL_t_mass_chi2: -3.15462
#gaussian_loglikelihood: obs=0.1181, mu=0.11772635, sigma=0.0011, sigma_th=0
#  Likelihood contribution from PrecisionBit::lnL_alpha_s_chi2: 5.83581

import numpy as np
import scipy.stats as sps
from functools import partial
import JMCtools.distributions as jtd
from JMCtools.experiment import Experiment

name  = ["top_mass", "alpha_s"]
obs   = [173.34, 0.1181]
sigma = [0.76,   0.0011]

def pars(loc,scale):
    return {"loc": loc, "scale": scale}

def get_seeds_full(samples,signal):
   loc = samples[...,0]
   return {'loc': loc} # We are directly sampling the MLEs, so this is trivial

def get_seeds_null(samples,signal):
   return {} # No nuisance parameters, so no nuisance parameter seeds

nuis_options = {} # None, no nuisance fit necessary

experiments = []
for n,o,s in zip(name,obs,sigma):
    joint = jtd.JointDist([jtd.TransDist(sps.norm,partial(pars,scale=s))])

    # Define the experiment object and options for fitting during statistical tests
    e = Experiment(n,joint,[o],DOF=1)

    general_options = {'loc': o, 'error_loc': s} # No real need for this either since seeds give exact MLE already.

    # For now we only define a 'gof' test, since there is no clear notion of a BSM contribution for these observables. At least not one that we can extract from our scan output.
    e.define_gof_test(null_options=nuis_options,
                  full_options=general_options,
                  null_seeds=(get_seeds_null, True), # extra flag indicates that seeds are exact
                  full_seeds=(get_seeds_full, True),
                  diagnostics=None
                  )

    experiments += [e]

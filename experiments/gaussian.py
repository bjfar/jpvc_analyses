"""Gaussian constraints that appear in MSSMEW.yaml"""

#gaussian_loglikelihood: obs=173.34, mu=175.04285, sigma=0.76, sigma_th=0
#  Likelihood contribution from PrecisionBit::lnL_t_mass_chi2: -3.15462
#gaussian_loglikelihood: obs=0.1181, mu=0.11772635, sigma=0.0011, sigma_th=0
#  Likelihood contribution from PrecisionBit::lnL_alpha_s_chi2: 5.83581
#  Likelihood contribution from ColliderBit::OPAL_Degenerate_Chargino_Conservative_LLike: 0
#gaussian_loglikelihood: obs=0.2, mu=0.19633835, sigma=0.0006, sigma_th=0.00010448813
#  Likelihood contribution from DecayBit::lnL_Z_invisible_width: -11.589

import JMCtools.distributions as jtd
import JMCtools.models as jtm
import numpy as np
import scipy.stats as sps
from functools import partial
from .experiment import Experiment

name  = ["top_mass", "alpha_s", "Z_invisible_width"]
obs   = [173.34, 0.1181, 0.2]
sigma = [0.76,   0.0011, 0.0006]

# s_MLE is not required, just useful for testing
s_MLE = obs # just cherry pick these for testing

def pars(loc,scale):
    return {"loc": loc, "scale": scale}

experiments = []
for n,o,s in zip(name,obs,sigma):
    e = Experiment(n)
    e.s_MLE = o #TESTING ONLY

    # Create the joint PDF object
    e.general_model = jtm.ParameterModel([jtd.TransDist(sps.norm,partial(pars, scale=s))])

    # Create the "observed" data
    # Need extra axes for matching shape of many simulated datasets
    e.observed_data = np.array([o])[np.newaxis,np.newaxis,:]
    
    # Define the null hypothesis
    # For general_model, this is the e.g. MSSM9 best fit value. We test this against the "free" alternate hypothesis.
    # For now this is just set such that the value tested exactly matches the observed value (so p=1, perfect fit)
    e.null_parameters = {'loc':o}
    
    # Define functions to get good starting guesses for fitting simulated data
    def get_seeds(samples):
        loc = samples[...,0]
        return {'loc': loc} # MLE's should just be the observed data
    
    # Same for fitting null hypothesis. But for this model there are no nuisance
    # parameters so there is nothing to fit
    def get_seeds_null(samples):
        return {}
    
    e.get_seeds = get_seeds
    e.get_seeds_null = get_seeds_null

    # Set options for fit
    # For null hypothesis there is nothing to fit! Let's see if we can force Minuit to
    # just compute the pdf for one point for us.
    e.null_options = {**e.null_parameters, 'fix_loc': True}
    e.general_options = {'loc': o, 'error_loc': s} #, 'print_level': 1} # Used observed value and std as starting guess and step size 
    e.nuis_options = {}
    
    # Degrees of freedom. Free parameters (not including mu) minus nuisance parameters
    e.DOF = 1
   
    # This won't make sense until we define a better null hypothesis, i.e. some SM-only hypothesis.
    # # Given a signal hypothesis, create a model we can test for effect size using just
    # # the free parameter 'mu'
    # def make_mu_model(s):
    #     # Create new parameter mapping functions with 'mu1' and 'mu2' parameters fixed.
    #     # The 'partial' tool from functools is super useful for this.
    #     s_model = jtm.ParameterModel([jtd.TransDist(sps.norm, partial(pars1, mu1=s[0])),
    #                                   jtd.TransDist(sps.norm, partial(pars2, mu2=s[1]))]
    #                                  ,[['mu'],['mu']])
    #     return s_model 

    def make_mu_model(s):
        pass
    e.make_mu_model = make_mu_model
   
    # Tell analysis code which sorts of models are defined for this experiment
    e.has_gof_model = True # Goodness of fit model ('general')
    e.has_mu_model = False  # 'mu' signal strength model ('mu_model')
 
    experiments += [e]

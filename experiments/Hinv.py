"""Higgs invisible width likelihood

   Here we essentially treat this the same as the simple Gaussian likelihoods,
   however it is a little more subtle because we need to infer the variance
   of the MLE from delta(chi^2) curves provided by the experiments
   and other fitting groups.

   This relies on the asymptotic normality of the MLE and the quadratic form
   of the log-likelihood surface.
"""

import JMCtools.distributions as jtd
import JMCtools.models as jtm
from JMCtools.experiment import Experiment
import scipy.stats as sps
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

name = "Higgs_invisible_width" 

# Data for Higgs invisible width likelihood
d1="DecayBit/data/arXiv_1306.2941_Figure_8.dat"
d2="DecayBit/data/CMS-PAS-HIG-17-023_Figure_7-b.dat"
chi2curve = np.loadtxt("/home/farmer/repos/gambit/copy3/{0}".format(d2))

# Suppose we only have this curve. How can we cook up a reasonable sampling
# distribution from this? I think there is little choice but to make some
# asymptotic assumptions about the distribution of the MLE (i.e. Gaussian)
# and fit to this curve to get the width.

# Suppose MLE is Gaussian in this parameter. The plot is for -2*delta(chi^2), by
# which they mean -2*log(L(BF)/L(\hat{BF})), presumably disallowing negative
# BF. In the asymptotic limit this curve is then
#
# (\hat{BF} - BF)^2 / sigma^2 + K
#
# So if we fit the curve to this, we can extract sigma, the standard deviation
# the MLE.
#
# Need function that computes residuals to do least squares fit
def chi2f(mu,muBF,sigma,K):
    return (muBF - mu)**2/sigma**2 + K

def res(x):
    muBF,sigma,K = x[0], x[1], x[2]
    x,y = chi2curve.T
    return chi2f(x,muBF,sigma,K) - y # Difference compared to data

x0 = np.array([0, 0.1, 0])
r = least_squares(res, x0)
hatBF, sigma, K = r.x

BF = np.arange(0,1,0.001)
chi2 = chi2f(BF, hatBF, sigma, K)
chi2_min = np.min(chi2) #Minimum over range [0,1]
dchi2 = chi2 - chi2_min

# Ok now build the probabilistic model for the MLE
def pars(BF):
    return {"loc":BF, "scale":sigma}

joint = jtd.JointDist([jtd.TransDist(sps.norm,pars)])

def get_seeds_full(samples,signal):
   BF = samples[...,0]
   return {'BF': BF} # We are directly sampling the MLEs, so this is trivial

def get_seeds_null(samples,signal):
   return {} # No nuisance parameters, so no nuisance parameter seeds

def get_asimov(mu,signal=None):
   # Need to return data for which mu=1 or mu=0 is the MLE
   BF = signal['BF']
   nA = mu*BF # I guess it is just this
   nuis_MLEs = {} # No nuisance parameters
   return nA, nuis_MLEs
 
nuis_options = {} # None, no nuisance fit necessary
general_options = {'BF': 0, 'error_BF': sigma} # No real need for this either since seeds give exact MLEs already.

#----------------------------------------------------------
# Setup done, now define the pdf and Experiment object

observed_data = np.array([hatBF]) 

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
                 scale_with_mu=['BF'],
                 )

e.define_musb_test(null_options=nuis_options,
                   mu1_seeds=(get_seeds_null,True),
                   mu0_seeds=(get_seeds_null,True),
                   scale_with_mu=['BF'],
                   asimov=get_asimov
                   )

experiments = [e]

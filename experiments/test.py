"""Simple Gaussian test experiment definitions

   We define a bunch of similar experiments, so
   that we can test their joint distributions.
"""

import JMCtools.distributions as jtd
import JMCtools.models as jtm
import numpy as np
import scipy.stats as sps
from functools import partial
from .experiment import Experiment

experiments = {}

# Parameter mapping function
def pars(s,b):
    return {"loc": s + b, "scale": 1}

def null_seeds(samples,signal):
    return {} # No nuisance parameters

def full_seeds(samples,signal,b):
    x = samples[:,0,0]
    #print("s:", x - b)
    return {"s": x - b} # Exact MLE for s

N = 10 # Number of Gaussian experiments to construct
for i in range(N):
        # We will set the "background" differently for each piece so we can tell them apart easier
        b = 20 + 5*i
        gauss = jtd.TransDist(sps.norm,partial(pars,b=b))

        # Create the "joint" PDF object (not very interesting since just one component)
        joint = jtd.JointDist([gauss])
         
        # Set options for parameter fitting
        s_opt  = {'s': 0, 'error_s': 1} # Will actually use seeds to obtain better starting guesses than this (actually, exact "guesses")
        
        nuis_options = {}  # No nuisance parameters (for now)
        general_options = {**s_opt}

        # Full observed data list, included observed values of nuisance measurements
        observed_data = [b + 5] # let's try a slight excess

        # Define the experiment object and options for fitting during statistical tests
        e = Experiment("gauss_{0}".format(i),joint,observed_data,DOF=1)
         
        e.define_gof_test(test_pars={'s': 0}, # Just for testing purposes
                          null_options=nuis_options,
                          full_options=general_options,
                          null_seeds=(null_seeds, True),
                          full_seeds=(partial(full_seeds,b=b), True),
                          )
        
        e.define_mu_test(nuisance_par_null={'s': 5}, # Cherry-picked to match observed data. Should give largest "local" p-value
                         null_options=nuis_options,
                         null_seeds=(null_seeds, True),
                         scale_with_mu=['s'], # 's' parameter to be scaled with signal strength 'mu' parameter
                         test_signal={'s': 5} # This should be overriden by some user choice, except for testing 
                         )
     
        experiments[e.name] = e

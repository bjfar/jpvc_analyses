"""Some toy modelling of the behaviour of look-elsewhere corrections

   I have some unanswered questions about when confidence regions can
   effectively tell us the level at which the null hypothesis can
   be rejected. So I want to compare these to explicit computations
   of global p-values.

   I have some intuition that in the large signal regime, the confidence
   intervals should give correct global p-values (if the parameter space
   is sufficiently regular in the vicinity of the best-fit region).

   This is based on the idea that asymptotic conditions should be better
   met if there is a large signal. In additon, the problem of parameters
   only existing under the alternate hypothesis (e.g. MSSM parameters
   don't exist for SM null hypothesis) is less of an issue (I think) if
   confidence intervals encapsulate a parameter region that is well-
   separated from the null-hypothesis-like boundary. 

   But that is all speculation. I want to test it numerically.
"""

import numpy as np
from JMCtools.analysis import Analysis
from JMCtools.experiment import Experiment

# We will use a simple model, something like the search for the SM Higgs boson.
# That is, let's have a Poisson likelihood, with some Gaussian width in a parameter
# "m", plus a strength parameter "mu", where the mass parameter "m" effectively
# vanishes for mu=0.

# ---- Create statistical model ----

def custpois(func,rename):
    """Construction transformed Poisson distribution,
       with logpmf function replaced by a version that
       can be evaluated for non-integer data (needed
       for Asimov calculationd"""

    mypois = jtd.TransDist(sps.poisson) # Null transformation, just to build object
    mypois.set_logpdf(smooth_poisson) # replace pdf calculation

    # Now build the transformed object
    return jtd.TransDist(mypois,func,rename)

def poisson_pars(s, b):
    l = np.atleast_1d(s + b)
    return {'mu': l}

Nbins = 100
background = 5 # Constant background in all bins
poisson_part = [custpois(partial(poisson_f_add,b=background),
                               ['s_{0} -> s'.format(i)], 
                         for i in range(Nbins)]

def get_seeds_null(samples,signal):
   """Seeds for nuisance parameter fits"""
   return {} # No nuisance parameters

def get_seeds(samples,signal):
   """Gets seeds for free bin signal parameter fits"""
   bin_samples = samples[:,0,:Nbins].T
   seeds = {}
   for i,x in enumerate(range(Nbins)):
      s_MLE = x - background
      seeds['s_{0}'.format(i)] = s_MLE
   return seeds

def get_asimov_data(mu,signal=None):
    if mu==0: # background-only case
        nA = np.ones(Nbins)*background
    else: # signal + background case; usually only want mu=1 but I allowed the more general case
        s = [signal['s_{0}'.format(i)] for i in range(self.N_SR)]
        nA = np.ones(Nbins)*background + mu*np.array(s)
    nuis_MLEs = {} # No nuisance parameters
    return nA, nuis_MLEs

# Create the joint PDF object
joint = jtd.JointDist(poisson_part)
         
s_opt  = {'s_{0}'.format(i): 0 for i in range(self.N_SR)} # Maybe zero is a good starting guess? Should use seeds that guess based on data.
        
nuis_options = {} # No nuisance pars
general_options = {**s_options, **nuis_options}

# Full observed data list, included observed values of nuisance measurements
observed_data = np.zeros(Nbins) # Will replace this with simulated data anyway

# Define the experiment object and options for fitting during statistical tests
e = Experiment(self.name,joint,observed_data,DOF=Nbins)
    
e.define_gof_test(null_options=nuis_options,
                  full_options=general_options,
                  null_seeds=(get_seeds_null, True),
                  full_seeds=(get_seeds, True), 
                  )

e.define_mu_test(null_options=nuis_options,
                 null_seeds=(get_seeds_null, False),
                 scale_with_mu=list(s_opt.keys()),
                 )

e.define_musb_test(null_options=nuis_options,
                   mu1_seeds=(self.seeds_null_f_gof(mu=1), False), # naming a bit odd, but these are the mu=1 seeds
                   mu0_seeds=(self.seeds_null_f_gof(mu=0), False), # " "   mu=0
                   scale_with_mu=list(s_opt.keys()),
                   asimov=get_asimov_data
                   )

tag = "gpval_properties"
expts = [e] # Only one experiment here
a = Analysis(expts,tag,make_plots=False)

# Set up signals to be tested

x = np.linspace(0,Nbins,1)


 



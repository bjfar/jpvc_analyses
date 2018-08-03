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
import scipy as sp
import scipy.stats as sps
import pickle
from functools import partial
import JMCtools.distributions as jtd
from JMCtools.analysis import Analysis
from JMCtools.experiment import Experiment
import matplotlib.pyplot as plt

# We will use a simple model, something like the search for the SM Higgs boson.
# That is, let's have a Poisson likelihood, with some Gaussian width in a parameter
# "m", plus a strength parameter "mu", where the mass parameter "m" effectively
# vanishes for mu=0.

# ---- Create statistical model ----

def smooth_poisson(k,mu,loc=0):
    r =  k*np.log(mu) - mu - sp.special.gammaln(k+1)
    return r

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
background = 20 # Constant background in all bins
poisson_part = [custpois(partial(poisson_pars,b=background),
                               ['s_{0} -> s'.format(i)])
                for i in range(Nbins)]

def get_seeds_null(samples,signal):
   """Seeds for nuisance parameter fits"""
   return {} # No nuisance parameters

def get_seeds(samples,signal):
   """Gets seeds for free bin signal parameter fits"""
   bin_samples = samples[:,0,:Nbins].T
   seeds = {}
   for i,x in enumerate(bin_samples):
      s_MLE = x - background
      seeds['s_{0}'.format(i)] = s_MLE
   #print("samples:", samples)
   #print("seeds:", seeds)
   return seeds

def get_asimov_data(mu,signal=None):
    if mu==0: # background-only case
        nA = np.ones(Nbins)*background
    else: # signal + background case; usually only want mu=1 but I allowed the more general case
        s = [signal['s_{0}'.format(i)] for i in range(Nbins)]
        nA = np.ones(Nbins)*background + mu*np.array(s)
    nuis_MLEs = {} # No nuisance parameters
    return nA, nuis_MLEs

# Create the joint PDF object
joint = jtd.JointDist(poisson_part)
         
s_opt  = {'s_{0}'.format(i): 0 for i in range(Nbins)} # Maybe zero is a good starting guess? Should use seeds that guess based on data.
s_opt2 = {'error_s_{0}'.format(i) :  np.sqrt(background) for i in range(Nbins)} # Get good step sizes from background fluctuation size
 
nuis_options = {} # No nuisance pars
general_options = {**s_opt, **s_opt2}

# Full observed data list, included observed values of nuisance measurements
observed_data = np.zeros(Nbins) # Will replace this with simulated data anyway

# Define the experiment object and options for fitting during statistical tests
name = "Toy_Higgs_search"
e = Experiment(name,joint,observed_data,DOF=Nbins)
    
e.define_gof_test(null_options=nuis_options,
                  full_options=general_options,
                  null_seeds=(get_seeds_null, True),
                  full_seeds=(get_seeds, True), 
                  )

e.define_mu_test(null_options=nuis_options,
                 null_seeds=(get_seeds_null, True),
                 scale_with_mu=list(s_opt.keys()),
                 )

e.define_musb_test(null_options=nuis_options,
                   mu1_seeds=(get_seeds_null, True), # naming a bit odd, but these are the mu=1 seeds
                   mu0_seeds=(get_seeds_null, True), # " "   mu=0
                   scale_with_mu=list(s_opt.keys()),
                   asimov=get_asimov_data
                   )

tag = "gpval_properties"
expts = [e] # Only one experiment here
a = Analysis(expts,tag,make_plots=False)

# Set up signals to be tested

# "energy range"
x = np.linspace(0,100,Nbins+1)

# Signal shape
def rate(bin_edges,cs,m):
    # Integrate by trapezoid method
    centres = 0.5*(bin_edges[1:] + bin_edges[0:-1])
    width = 5
    s = cs*sps.norm.pdf(x=centres,loc=m,scale=width) # times some bin width that doesn't matter
    return s

# Get grid of signal predictions over (cs,m) plane

cs = np.linspace(0,200,20)
m = np.linspace(0,100,40)

r = rate(x,cs[:,np.newaxis,np.newaxis],m[:,np.newaxis])
CS, M = np.meshgrid(cs,m,indexing='ij')

print(r.shape)
print(CS.shape)
print(M.shape)

# # Plot rates to check we did it right
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(x[:-1],r[25,50,:]+background,drawstyle='steps-post',label="cs={0:.2f}, m={1:.2f}".format(cs[25],m[50]))
# ax.plot(x[:-1],r[10,50,:]+background,drawstyle='steps-post',label="cs={0:.2f}, m={1:.2f}".format(cs[10],m[50]))
# ax.plot(x[:-1],r[10,25,:]+background,drawstyle='steps-post',label="cs={0:.2f}, m={1:.2f}".format(cs[10],m[25]))
# ax.set_xlabel("Energy bins")
# ax.set_ylabel("signal rate")
# ax.legend(loc=1, frameon=False, framealpha=0,prop={'size':10})
# fig.savefig("bump_hunt_test.png")

#CS, M = np.mgrid[0:50:1j, 0:100:1j]
#positions = np.vstack([X.ravel(), Y.ravel()])

# Convert rate array into list of signal parameters
signals = []
for point in r.reshape(-1,r.shape[-1]):
    signals += [{"s_{0}".format(i):s for i,s in enumerate(point)}]

# Generate some data to use as the observed data
Nsamples = 1e4
test_type = 'musb'
#true_r = rate(x,cs=70,m=60) # Fake signal to inject
true_r = rate(x,cs=0,m=0) # Fake signal to inject
true_parameters = {name: {"s_{0}".format(i):s for i,s in enumerate(true_r)}}
obs = a.simulate(Nsamples,test_type,true_parameters)[0]

# Evaluate likelihood function for these signal hypotheses
# This actually gets -2*log(L_s+b/L_b), but L_b is common to all so is a constant factor anyway.
N = len(signals)
LLRs = np.zeros((N,Nsamples))
apvals = np.ones((N,Nsamples)) # Easiest to save them all, and THEN calculate the minima over parameter points
gof_LLRs = np.zeros((N,Nsamples))
gof_apval = np.zeros((N,Nsamples))

# Do vanilla GOF test for background-only hypothesis
model, gof_LLR, gof_LLR_obs, gof_apval, gof_epval, gofDOF = e.do_gof_test(s_opt,samples=None,observed=obs)

# Do musb test for all signal hypotheses
for i,sig in enumerate(signals):
    print("Getting LLR for point {0}".format(i))
    model, musb_LLR, musb_LLR_obs, musb_apval, musb_epval, LLRA, Eq, Varq = e.do_musb_test(sig,samples=None,nullmu=0,observed=obs)
    LLRs[i] = musb_LLR_obs
    apvals[i] = musb_apval

# pickle results
with open("toy_higgs_{0}.pkl".format(tag), 'wb') as pkl_file: 
    pickle.dump([background,x,obs,CS,M,r,LLRs,apvals,gof_LLR_obs,gof_apval,gofDOF],pkl_file)



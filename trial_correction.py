"""Trial-factor corrections for MSSMEW4 scan

   For now; asymptotic p-values only (but these
   are pretty accurate it seems)

   For this computation, we take pseudo-data
   generated under the background-only hypothesis,
   and, under each pseudo-dataset, compute
   "observed" local p-values for all sampled points,
   saving the minimum p-value found.
   We then infer a global p-value based on the
   simulated distribution of minimum local p-values.
"""

import numpy as np
import scipy.stats as sps
import h5py
import six
from JMCtools.analysis import Analysis
from JMCtools.experiment import Experiment
import experiments.CBit_LHC_python as CBa
import datatools as dt

# Input nominal signal predictions for the
# Collider experiments to be analysed

# Will extract a whole bunch of these from MSSMEW run
#hdf5file = "/home/farmer/repos/gambit/copy3/runs/MSSMEW/samples/MSSMEW.hdf5"
hdf5file = "data/MSSMEW_pp_bestfit_region.hdf5"
gname = "MSSMEW"
f = h5py.File(hdf5file,'r')
g = f[gname]

# Signal predictions in hdf5 file have dataset names like this: 
#LHC_signals\ @ColliderBit::calc_LHC_signals::ATLAS_13TeV_3b_24invfb__meff340_ETmiss70__i14__signal
#LHC_signals\ @ColliderBit::calc_LHC_signals::ATLAS_13TeV_3b_24invfb__meff340_ETmiss70__i14__signal_uncert

# Need to load these predictions up and assign them to the corresponding Python analyses.
# The names are at least systematic so we can avoid repeating most of the junk.
# We name the Python analyses to match the GAMBIT ones, as well as the signal region names.

hdf5_names = {} # Names of hdf5 datasets to load
for a in CBa.analyses:
   hdf5_names[a.name] = ["#LHC_signals @ColliderBit::calc_LHC_signals::{0}__{1}__signal".format(a.name,SR) for SR in a.SR_names]

def get_signal(aname,m,i):
    """Extract signal vector for analysis 'name', chosen by applying mask 'm' 
       and then index 'i' to the corresponding HDF5 datasets"""
    dsets = dt.get_data(g, hdf5_names[aname], m, i)
    return [d.data() for d in dsets] 

tag = "GlobalTest_1e3"
Nsamples = int(1e3)

# Choose which points to analyses
logl = dt.get_data(g, ["LogLike"])[0]
m = logl.valid() # all points with valid likelihoods
#N = np.sum(m)
N = 10 # testing

# Begin loop over signal hypotheses in HDF5 file
# 
pvals = np.ones((N,Nsamples)) # Easiest to save them all, and THEN calculate the minima over parameter points
pseudodata_b = None # Only generate the pseudodata once; use same samples for all parameter points
for i in range(N):
    if N==1: i = None
    # Extract the signal predictions for each point
    print("Analysing HDF5 entry with signal:")
    CBsignal = {}
    for a in CBa.analyses:
        CBsignal[a.name] = get_signal(a.name,m,i)
        print("  {0}: {1}".format(a.name, CBsignal[a.name]))

    # Generate experiment objects to analyse
    expts = []
    test_parameters = {} # Signal predictions of point in parameter space to be tested
    true_gof_parameters = {} # Signal + nuisance parameters of point in parameters space used to generate data (always defined in "gof" parameter space)
    true_musb_parameters = {}
    
    for a in CBa.analyses:
        # Signal hypothesis needs to be supplied prior to building Experiment 
        # objects so that we can use it to select which signal regions to use for 
        # the analysis (as in ColliderBit)
        s = CBsignal[a.name]
        s_dict = {'s_{0}'.format(i): val for i,val in enumerate(s)}
        b_dict = {'s_{0}'.format(i): 0 for i in range(len(s))}
        nuis_dict = {'theta_{0}'.format(i): 0 for i in range(len(s))} # nominal values of nuisance parameters (for simulating data)
        test_parameters[a.name] = s_dict
        true_gof_parameters[a.name]  = {**s_dict, **nuis_dict}
        true_musb_parameters[a.name] = {**b_dict, **nuis_dict}
        try:
           expts += [a.make_experiment(s_dict)]
        except:
           six.reraise(Exception,Exception("Mystery error encountered while calling make_experiment for analysis {0}".format(a.name))) 

    cb = Analysis(expts,tag)
   
    # Generate pseudodata
    if pseudodata_b is None:
       pseudodata_b = cb.simulate(Nsamples,'musb',true_musb_parameters) #background-only pseudodata this time)

    # Test significance of data under background-only hypothesis (with signal at "test_parameters" being the signal hypothesis)
    cb.musb_analysis(test_parameters,pseudodata=None,nullmu=0,observed=pseudodata_b) # pseudodata=None means use asymptotic theory to compute p-values
 
     
    mr = cb.results("(test == 'musb_mu=0') and (experiment == 'Monster')")
    print(mr)
    p = mr["a_pval"].item()
    print(p)
    pvals[i] = p 
    # how to save results? Would be nice to see component-wise, i.e. trial-corrected values for all analyses as well as the combination

print("All p-values:")
print(p)

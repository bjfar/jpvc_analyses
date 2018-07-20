"""p-value calculations for MSSMEW4 scan"""

import numpy as np
import scipy.stats as sps
import pandas as pd
import h5py
import six
import pickle
from JMCtools.analysis import Analysis
from JMCtools.experiment import Experiment
import experiments.CBit_LHC_python as CBa # ColliderBit analyses
import experiments.Hinv as Hinv # Higgs invisible width
import experiments.Zinv as Zinv # Z invisible width
import experiments.gaussian as gaussian # alpha_s and m_top
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

tag = "AllLikes_1e4"
Nsamples = int(1e4)

# Choose which points to analyses
logl = dt.get_data(g, ["LogLike"])[0]
logl_max = np.max(logl.valid())
m = (logl == logl_max) # Selecting just the maximum likelihood point
N = 1

# Skip the analysis and just show us the signal predictions to be analysed
signal_only = False

# Begin loop over signal hypotheses in HDF5 file
for i in range(N):
    if N==1: i = None

    # Generate storage for experiment objects to be analysed,
    # and parameters to use for various tests
    expts = []
    test_parameters = {} # Signal predictions of point in parameter space to be tested
    mu1_parameters = {} # Signal + nuisance parameters of point in parameters space used to generate data (always defined in "gof" parameter space)
    mu0_parameters = {}
    gof_test_parameters = {} # Signal parameters to use for gof test hypothesis (nuisance left free to float)
    gof_parameters = {} # Parameters used to generate pseudodata for gof test
    gof_b_test_parameters = {}
    gof_b_parameters = {}     


    # Collect experiments to be analysed
 
    # ------ Higgs invisible width:
    expts += Hinv.experiments
    dset = "#inv_Higgs_BF @DecayBit::MSSM_inv_Higgs_BF"
    BF = dt.get_data(g, [dset], m, i)[0].data()
    test_parameters    ["Higgs_invisible_width"] = {'BF': BF}  
    gof_test_parameters["Higgs_invisible_width"] = {'BF': BF}  
    gof_b_test_parameters["Higgs_invisible_width"] = {'BF': 0}  
    mu1_parameters     ["Higgs_invisible_width"] = {'BF': BF}  
    mu0_parameters     ["Higgs_invisible_width"] = {'BF': 0}  
    gof_parameters     ["Higgs_invisible_width"] = {'BF': BF}  
    gof_b_parameters   ["Higgs_invisible_width"] = {'BF': 0}  
    # ======

    # ------ Z invisible width:
    expts += Zinv.experiments
    # SM prediction varies with nuisance parameters, so need to extract this 
    # point by point. I guess ideally we'd profile over those too, but we
    # can't do that because everything is pre-computed.
    dsets = [
"#Z_gamma_chi_0 @DecayBit::Z_gamma_chi_0_MSSM_tree::central",
"#Z_gamma_chi_0 @DecayBit::Z_gamma_chi_0_MSSM_tree::lower",
"#Z_gamma_chi_0 @DecayBit::Z_gamma_chi_0_MSSM_tree::upper",
"#Z_gamma_nu @DecayBit::Z_gamma_nu_SM_2l::central",
"#Z_gamma_nu @DecayBit::Z_gamma_nu_SM_2l::lower",
"#Z_gamma_nu @DecayBit::Z_gamma_nu_SM_2l::upper",
]
    gamma_inv_BSM, gamma_inv_BSM_upper, gamma_inv_BSM_lower,\
    gamma_inv_SM, gamma_inv_SM_upper, gamma_inv_SM_lower\
      = dt.get_data(g, dsets, m, i)

    print("gamma_inv_BSM:", gamma_inv_BSM.data())
    # Ok for my test point this is zero, so musb can do nothing.
    # For testing we will just set this to something non-zero
    gamma_inv_BSM = gamma_inv_BSM.data()
    #gamma_inv_BSM = -0.001 # The SM contribution is actually slightly larger than measured, so we need to suppress it to get a better fit

    # Average + and - errors (this is what is done in GAMBIT)
    tau_SM = 0.5 * (gamma_inv_SM_upper.data() + gamma_inv_SM_lower.data());
    tau_BSM = 0.5 * (gamma_inv_BSM_upper.data() + gamma_inv_BSM_lower.data());
    # Add theory errors in quadrature (again as done in GAMBIT)
    tau = np.sqrt(tau_SM**2 + tau_BSM**2)
    common_pars = {"gamma_inv_SM": gamma_inv_SM.data(),
                   "sigma_err": tau}
    test_parameters    ["Z_invisible_width"] = {'gamma_inv_BSM': gamma_inv_BSM, **common_pars}  
    gof_test_parameters["Z_invisible_width"] = {'gamma_inv_BSM': gamma_inv_BSM, **common_pars}  
    gof_b_test_parameters["Z_invisible_width"]={'gamma_inv_BSM': 0, **common_pars}  
    mu1_parameters     ["Z_invisible_width"] = {'gamma_inv_BSM': gamma_inv_BSM, **common_pars}  
    mu0_parameters     ["Z_invisible_width"] = {'gamma_inv_BSM': 0, **common_pars}  
    gof_parameters     ["Z_invisible_width"] = {'gamma_inv_BSM': gamma_inv_BSM, **common_pars}  
    gof_b_parameters   ["Z_invisible_width"] = {'gamma_inv_BSM': 0, **common_pars}  
    # ======

    # ------ Simple Gaussian constraints (alpha_s and top mass)
    expts += gaussian.experiments 
    dsets = [
"#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::alphaS",
"#StandardModel_SLHA2_parameters @StandardModel_SLHA2::primary_parameters::mT"
]
    alpha_S, mtop = dt.get_data(g, dsets, m, i)

    gof_test_parameters["alpha_s"]  = {'loc': alpha_S.data()}  
    gof_test_parameters["top_mass"] = {'loc': mtop.data()}  
    gof_parameters["alpha_s"]  = {'loc': alpha_S.data()}  
    gof_parameters["top_mass"] = {'loc': mtop.data()}  
    gof_b_test_parameters["alpha_s"]  = {'loc': alpha_S.data()}  
    gof_b_test_parameters["top_mass"] = {'loc': mtop.data()}  
    gof_b_parameters["alpha_s"]  = {'loc': alpha_S.data()}  
    gof_b_parameters["top_mass"] = {'loc': mtop.data()}  
     # ======  

    # ------ ColliderBit experiments:
    # Extract the signal predictions for each point
    print("Analysing HDF5 entry with signal:")
    CBsignal = {}
    for a in CBa.analyses:
        CBsignal[a.name] = get_signal(a.name,m,i)
        print("  {0}: {1}".format(a.name, CBsignal[a.name]))
    if signal_only: 
        continue
   
    for a in CBa.analyses:
        # Signal hypothesis needs to be supplied prior to building Experiment 
        # objects so that we can use it to select which signal regions to use for 
        # the analysis (as in ColliderBit)
        s = CBsignal[a.name]
        s_dict = {'s_{0}'.format(i): val for i,val in enumerate(s)}
        b_dict = {'s_{0}'.format(i): 0 for i in range(len(s))}
        nuis_dict = {'theta_{0}'.format(i): 0 for i in range(len(s))} # nominal values of nuisance parameters (for simulating data)
        test_parameters[a.name] = s_dict
        gof_test_parameters[a.name] = s_dict # Testing goodness of fit of signal+background hypothesis 
        gof_b_test_parameters[a.name] = b_dict # Testing goodness of fit of background-only hypothesis 
        mu1_parameters[a.name]  = {**s_dict, **nuis_dict}
        mu0_parameters[a.name]  = {**b_dict, **nuis_dict}
        gof_parameters[a.name]  = {**s_dict, **nuis_dict}
        gof_b_parameters[a.name]  = {**b_dict, **nuis_dict}
        try:
           expts += [a.make_experiment(s_dict)]
        except:
           six.reraise(Exception,Exception("Mystery error encountered while calling make_experiment for analysis {0}".format(a.name))) 
    # ======


    # ----- Run analysis

    cb = Analysis(expts,tag)
    # Need pseudodata under both mu=1 and mu=0 hypotheses
    pseudodata_mu1 = cb.simulate(Nsamples,'musb',mu1_parameters)
    pseudodata_mu0 = cb.simulate(Nsamples,'musb',mu0_parameters)

    # Do goodness of fit test (signal at "test_parameters" vs best fit in general signal space)
    # Can also test goodness of fit of background-only model against fully-general signal
    #print("gof_test_parameters:", gof_test_parameters)
    pseudodata_gof_b = cb.simulate(Nsamples,'gof',gof_b_parameters)
    cb.gof_analysis(gof_b_test_parameters,pseudodata_gof_b)
    #print(cb.results())
    #quit()

    # Want to do MSSMEW4 goodness of fit too, so need to rename previous 'gof' results in outputpiif allow_corr:
    df = cb.results().results_df # Underlying Pandas dataframe
    df['test'].replace('gof','gof_b',inplace=True)

    # Now run 'gof' analysis again
    # Kind of wasteful to do all the fitting again, but currently I have no to mechanism to store the previous pdfs etc.
    pseudodata_gof = cb.simulate(Nsamples,'gof',gof_parameters)
    cb.gof_analysis(gof_test_parameters,pseudodata_gof)      

    # Test significance of data under signal hypothesis (with signal at "test_parameters" being the signal hypothesis)
    # Same pseudodata as before (signal hypothesis pseudodata)
    #cb.musb_analysis(test_parameters,pseudodata_mu1,nullmu=1)
     
    # Test significance of data under background-only hypothesis (with signal at "test_parameters" being the signal hypothesis)
    #cb.musb_analysis(test_parameters,pseudodata_mu0,nullmu=0)
    
    # Test both mu=0 and mu=1 at once
    cb.musb_analysis_dual(test_parameters,pseudodata_mu1,pseudodata_mu0)
 
    # Compute one-tailed significances from p-values
    ap = cb.results()["a_pval"]
    ep = cb.results()["e_pval"]
    #print("apvals:", ap)
    #print("epvals:", ep)
    asig = -sps.norm.ppf(ap)
    esig = -sps.norm.ppf(ep)
    cb.results().add_column("a_sig",asig)
    cb.results().add_column("e_sig",esig)

    # Present the results computed so far as a table
    # The object remembers what tests have been done, and stores the summary results internally.
    print(cb.results('test == "gof"')) 
    print(cb.results('test == "gof_b"')) 
    print(cb.results('test == "musb_mu=1"'))
    print(cb.results('test == "musb_mu=0"'))
    #print(cb.results('test == "musb_CLs"'))

    # Pickle results so we can more freely fiddle about with format for publication
    with open("MSSMEW4_pvalue_results.pkl", 'wb') as pkl_file: 
        pickle.dump(cb.results(),pkl_file)    


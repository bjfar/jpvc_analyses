from JMCtools.analysis import Analysis
from JMCtools.experiment import Experiment

import experiments.ColliderBit_analysis as CBa

# Input nominal signal predictions for the
# Collider experiments to be analysed

CBsignal = {} 
CBsignal["ATLAS_13TeV_MultiLEP_36invfb"]  = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0, 1, 4, ]
CBsignal["CMS_13TeV_1LEPbb_36invfb"]      = [0, 8, ]
CBsignal["CMS_13TeV_2LEPsoft_36invfb"]    = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, ]
# ###CBsignal["CMS_13TeV_2OSLEP_36invfb"]    = [0, 1, 3, 2, 0, 0, 0, ]
# ###CBsignal["CMS_13TeV_2OSLEP_confnote_36invfb_NOCOVAR_NOLIKE"]     = [0, 0, 0, 1, 1, 0, 0, 0, 0, ]
CBsignal["CMS_13TeV_MONOJET_36invfb"]     = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]
CBsignal["CMS_13TeV_MultiLEP_36invfb"]    = [0, 2, 1, 2, 1, 1, 0, ]

# Compile the list of experiments to analyse
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
    expts += [a.make_experiment(s_dict)]   
 
tag = "5e3"
Nsamples = int(float(tag))
 
cb = Analysis(expts,tag)
# Do goodness of fit test (signal at "test_parameters" vs best fit in general signal space)
pseudodata_gof = cb.simulate(Nsamples,'gof',true_gof_parameters)
#pseudodata_gof = None # testing asymptotic only
cb.gof_analysis(test_parameters,pseudodata_gof)

# Test significance of data under background-only hypothesis (with signal at "test_parameters" being the alternate hypothesis)
pseudodata_musb = cb.simulate(Nsamples,'musb',true_musb_parameters)
#pseudodata_musb = None # testing asymptotic only
cb.musb_analysis(test_parameters,pseudodata_musb)

# Present the results computed so far as a table
# The object remembers what tests have been done, and stores the summary results internally.
print(cb.results('test == "gof"'))

print(cb.results('test == "musb"'))

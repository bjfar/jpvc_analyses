import numpy as np
import scipy.stats as sps
import h5py
import pickle
import datatools as dt
import experiments.CBit_LHC_python as LHC
import CMS_likelihood_recompute as CMSlike
import plottools as t
import matplotlib.pyplot as plt

dfile = "data/MSSMEW_pp_final.hdf5"

gname = "MSSMEW"
f = h5py.File(dfile,'r')
g = f[gname]

# Analyses whose likelihoods we need to recompute
recompute_for = [
 "CMS_13TeV_2LEPsoft_36invfb"
#,"CMS_13TeV_2OSLEP_36invfb"
]

# Combined likelihoods (will need to re-weight these)
loglikes = {}
analyses = [a for a in LHC.analyses if a.name in recompute_for]
for a in analyses:
    loglikes[a.name] = "#LHC_LogLike_per_SR @ColliderBit::get_LHC_LogLike_per_SR::{0}__combined_LogLike".format(a.name)
LHC_loglike = "#LHC_Combined_LogLike @ColliderBit::calc_combined_LHC_LogLike"
LogLike = "LogLike"

llike = dt.get_data(g, [LogLike])[0] #, m, i)

N = len(llike.data())
#N = 100 # testing

# First we need the marginal SM likelihoods, i.e. with signal=0
logl_b = np.zeros((1,2))
X = {}
for i,a in enumerate(analyses):
    logl_b[0,i] = CMSlike.compute_marg_logl(a.SR_n,a.SR_b,a.cov,M=int(1e5))

print("logl_b:",logl_b)
quit()

# Now load up the MSSM4 likelihoods we computed earlier
with open("data/new_CMS_likes_finished.pkl", 'rb') as pkl_file: 
    CMSlikes = pickle.load(pkl_file)

print(CMSlikes.shape)
#for i in range(1000):
#    print(CMSlikes[i])

print(np.max(CMSlikes - logl_b,axis=0))

CMSlogl = CMSlikes - logl_b
 
# Let's make a few plots to check that these results make sense

makeplots =False
if makeplots:
    dsets = [
     "#MSSM_spectrum @SpecBit::get_MSSM_spectrum_as_map::~chi0_1 Pole_Mass"
    ,"#MSSM_spectrum @SpecBit::get_MSSM_spectrum_as_map::~chi+_1 Pole_Mass"
    ]
    
    mchi0, mchip = dt.get_data(g, dsets)
    print(mchi0.data().shape)
    print(mchip.data().shape)
    print(CMSlogl.shape)
    
    m = (np.abs(mchi0.data()) < 700) & (np.abs(mchip.data()) < 700)
    
    for i,a in enumerate(analyses):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        data = np.vstack([np.abs(mchip.data())[m],np.abs(mchi0.data())[m],-2*CMSlogl[:,i][m]]).T
        t.profplot(ax,data,title=a.name,labels=["mchi+","mchi0"],nybins=100,nxbins=100)
        fig.savefig("{0}_newlike_prof.png".format(a.name)) 
        plt.close(fig)
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        t.chi2scatplot(ax,data,title=a.name,labels=["mchi+","mchi0"])
        fig.savefig("{0}_newlike_scat.png".format(a.name)) 
        plt.close(fig)

# Cool, looks reasonable.
# Now we need to re-weight the various likelihoods which these results affect
# Which I think are these:
# loglikes[a.name] = "#LHC_LogLike_per_SR @ColliderBit::get_LHC_LogLike_per_SR::{0}__combined_LogLike".format(a.name)
# LHC_loglike = "#LHC_Combined_LogLike @ColliderBit::calc_combined_LHC_LogLike"
# LogLike = "LogLike"

# Open copy of results file for output
dfile_out = "outdata/MSSMEW_pp_final_reprocess.hdf5"
f_out = h5py.File(dfile_out,'r+')
g_out = f_out[gname]

#  old_logl_piece = np.zeros(llike.data().shape)
#  new_logl_piece = np.zeros(llike.data().shape)
for i,a in enumerate(analyses):
    #olddata = g[loglikes[a.name]]
    newdata = g_out[loglikes[a.name]]
    newdata[...] = CMSlogl[:,i]
    
    #old_logl_piece += olddata[:]
    #new_logl_piece += CMSlogl[:,i]
  
#  print("Old and new logl pieces")
#  print(np.sort(old_logl_piece)[::-1])
#  print(np.sort(new_logl_piece)[::-1])
#  
#  LHC_loglike = "#LHC_Combined_LogLike @ColliderBit::calc_combined_LHC_LogLike"
#  LogLike = "LogLike"
#  
#  LHCll_in  = g[LHC_loglike]
#  lltot_in  = g[LogLike]  
#  LHCll_out = g_out[LHC_loglike]
#  lltot_out = g_out[LogLike]
#  
#  LHCll_out[...] = LHCll_in - old_logl_piece# + new_logl_piece
#  lltot_out[...] = lltot_in - old_logl_piece# + new_logl_piece
#  #LHCll_out[...] = LHCll_in + new_logl_piece
#  #lltot_out[...] = lltot_in + new_logl_piece
#  
#  print("Old and new total logls")
#  print(np.sort(LHCll_in) [::-1])
#  print(np.sort(LHCll_out)[::-1])
#  print(np.sort(lltot_in) [::-1])
#  print(np.sort(lltot_out)[::-1])

# Copy the old loglike numbers to new datasets, so Anders can use them during reweighting
for i,a in enumerate(analyses):
    dset_name = loglikes[a.name]
    old_loglike = g[dset_name][:]
    old_loglike_valid = g[dset_name+"_isvalid"][:]
    g_out.create_dataset("{0}_old".format(dset_name), data=old_loglike, chunks=(1000,))
    g_out.create_dataset("{0}_old_isvalid".format(dset_name), data=old_loglike_valid, chunks=(1000,))


print("Done!")

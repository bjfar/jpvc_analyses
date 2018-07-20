"""Analysis of trial_correction.py output"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as spi
import scipy.stats as sps
import JMCtools.common as c

tag = "GlobalTest_1e4"
with open("LEEpvals_{0}.pkl".format(tag), 'rb') as pkl_file: 
    pvals = pickle.load(pkl_file)

Nsamples = pvals.shape[1]

print("All p-values:")
# dims are (N parameter points, M pseudodata samples)
# So we need to minimise over the pseudodata samples)
valid_pvals = pvals[np.all(pvals!=1.,axis=1)] # Generation may have been cut off; array is padded with ones.
print(valid_pvals.shape)

#print("Minimum p-values over all parameter points")
#print(np.min(pvals,axis=0))
minpvals = np.min(valid_pvals,axis=0)

print("len(pvals[0]):",len(pvals[0]))
print("len(minpvals):",len(minpvals))

minp = np.min(minpvals)
print("Smallest p-value observed:", minp)

print("Generating CDF plots")
minp = np.log10(minp) # 10^-4
x = np.logspace(minp,0,10)
# Ok let's make some plots of the CDF, to see the look-elsewhere effect in action
fig = plt.figure(figsize=(10,8))
nbins = int(0.1*10**(-minp))

# First, the local version
ax = fig.add_subplot(221)
#n, bins = np.histogram(valid_pvals[0], bins=nbins) #, normed=True)
#n = np.cumsum(n)/Nsamples
# We can do better than a histogram
s = np.sort(valid_pvals[0])
CDF = c.eCDF(s)
ax.plot(s,CDF,drawstyle='steps-post',label="Simulated",c='b')
ax.plot(x,x)
ax.set_xlabel("local p-value")
ax.set_ylabel("local cdf(local p-value)")
ax.set_xscale("log")     
ax.set_yscale("log")     
#ax.set_xlim(ran[0],ran[1])
#ax.set_ylim(yran[0],yran[1])
  
# Now the global version
ax = fig.add_subplot(222)
trials = Nsamples
#n, bins = np.histogram(minpvals, bins=nbins) #, normed=True)
#n = np.cumsum(n)/Nsamples
s = np.sort(minpvals)
CDF = c.eCDF(s)
ax.plot(s,CDF,drawstyle='steps-post',label="Simulated",c='b')
ax.plot(x,x)
ax.set_xlabel("min(local p-value)")
ax.set_ylabel("cdf(min(local p-value))")
ax.set_xscale("log")     
ax.set_yscale("log")     

# Global version with significance axes
ax = fig.add_subplot(223)
lsig = -sps.norm.ppf(s)
gsig = -sps.norm.ppf(CDF)
m = (lsig>0) & (gsig>0)
ax.plot(lsig[m],gsig[m],drawstyle='steps-post',label="Simulated",c='b')
y = np.linspace(0,5,10) 
ax.plot(y,y)
ax.set_xlabel("local significance")
ax.set_ylabel("global significance")

# Plot showing evolution of global significance as more parameter points are added (randomised selection)
ax = fig.add_subplot(224)
np.random.shuffle(valid_pvals)
curves = [[] for i in range(5)]
for i in range(1,len(valid_pvals)):
    minpvals = np.min(valid_pvals[:i],axis=0)
    s = np.sort(minpvals)
    CDF = c.eCDF(s)
    lsig = -sps.norm.ppf(s)
    gsig = -sps.norm.ppf(CDF)
    a = np.argsort(lsig)
    # Interpolate
    GSIG = spi.interp1d([0]+list(lsig[a])+[100],[0]+list(gsig[a])+[100])
    for i in range(5):
        curves[i] += [GSIG(i+1)] 
for i in range(5):
    ax.plot(range(1,len(valid_pvals)),curves[i],drawstyle='steps-post',label="local sigma = {0}".format(i+1))
ax.legend(loc=1, frameon=False, framealpha=0,prop={'size':10})
ax.set_xlabel("Number of parameter points tested")
ax.set_ylabel("global significance")

plt.tight_layout()
fig.savefig("LEE_cdf.png")



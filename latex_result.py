"""Converting results into nice LaTeX tables"""

import pandas as pd
import pickle
import numpy as np
import scipy.stats as sps
from functools import reduce

pd.set_option('display.max_colwidth', -1)

path = "outdata"
#tag = "AllLikes_BF"
#tag = "CMSfixed_AllLikes_BF"
#tag = "rev1_with_8TeV"
tag = "rev1_without_8TeV"

# Normal results, using best expected signal region
with open("{0}/MSSMEW4_pvalue_results_{1}_2e4.pkl".format(path,tag), 'rb') as pkl_file: 
    results = pickle.load(pkl_file)  

# Also get results using all signal regions, ignoring correlations when no covariance matrix is available
with open("{0}/MSSMEW4_pvalue_results_{1}_2e4_NOCORR.pkl".format(path,tag), 'rb') as pkl_file: 
    results_NOCORR = pickle.load(pkl_file)  

# # Extra no-correlation results for GOF test, with discovery SRs used for ATLAS 3b analysis
with open("{0}/MSSMEW4_pvalue_results_{1}_2e4_NOCORR_discoverySR.pkl".format(path,tag), 'rb') as pkl_file: 
    results_NOCORR_DSR = pickle.load(pkl_file)  

#print(results_NOCORR)
#quit()

# Prepare latex tables for publication
# Need to rearrange them a bit, and improve labels
# Want to merge the background-only gof and musb_mu=0 results
# So pull out the underlying pandas dataframes, then add new columns to one of them
s_gof = results.query('test == "gof"').results_df
bg_gof = results.query('test == "gof_b"').results_df
musb_0 = results.query('test == "musb_mu=0"').results_df
#s_gof_NC = results_NOCORR.query('test == "gof"').results_df # Gave up on using discovery signal regions for these, don't have them for the new 8 TeV analyses
#bg_gof_NC = results_NOCORR.query('test == "gof_b"').results_df
s_gof_NC = results_NOCORR_DSR.query('test == "gof"').results_df
bg_gof_NC = results_NOCORR_DSR.query('test == "gof_b"').results_df
musb_0_NC = results_NOCORR.query('test == "musb_mu=0"').results_df
dframes = [s_gof,bg_gof,musb_0,s_gof_NC,bg_gof_NC,musb_0_NC]

# Set experiment names as index
ename='Analysis'
for d in dframes:
    d.rename(columns={'experiment': ename}, inplace=True)
    d.set_index(ename, inplace=True)

# In the DSR datasets, rename the 3d "discovery" analysis so it merges with the 
# others. We will explain the difference in a footnote or something.
for d in [s_gof_NC,bg_gof_NC]:
    d.rename(index={"ATLAS_13TeV_3b_discoverySR_24invfb": "ATLAS_13TeV_3b_24invfb"}, inplace=True)

# --- NEW! Meta-analysis p-value combination ---
# Compute combined p-values obtained via a meta-analysis method
# (Fisher's method in this case)
# I am interested to see if this suffers less dilution due to
# bazillions of degrees of freedom in some signal regions

new_dframes = []
for df in dframes:
    test = df['test'].iloc[0]
    new_record = [test]
    DOF = None
    pvals = []
    sigs = []
    for pname in ['a_pval','e_pval']:
        p = df[pname].iloc[:-1] # leave the Monster pvalue out!!!
        # Compute combination
        x = -2*np.sum(np.log(p)) # Fishers method
        DOF = 2*len(p)
        p_comb = 1 - sps.chi2.cdf(x,df=DOF)
        sig_comb = -sps.norm.ppf(p_comb)
        #print("p_comb:", p_comb, "sig_comb:", sig_comb)
        pvals += [p_comb]
        sigs += [sig_comb]
    new_record += pvals + [DOF] + sigs
    # Add new record
    df.loc['Meta'] = new_record
    #pd.DataFrame([new_record],columns=list(df))
    #new_dframes += [pd.concat([df,new_df])]
    #["experiment","test","a_pval","e_pval","DOF","a_sig","e_sig"]

# ---

# Pick whether to use empirical or asymptotic results
use = 'e_sig' 
nouse = 'a_sig'

# # Delete all the stuff we don't want in the result
# EDIT: No need, just select what we want at the end
# del_list = ['test', 'a_pval', 'e_pval', nouse] # keep use
# for df in [bg_gof,musb_0]:
#     for col in del_list:
#         df.pop(col)

# Do some renaming
use_shortnames = True
if not use_shortnames:
    lname=r'\makecell{Local \\ significance ($\sigma$)}'
    lNCname=r'\makecell{Local \\ significance ($\sigma$) \\ w. no corr.}'
    gbname=r'\makecell{Background \\ goodness-of-fit ($\sigma$)}'
    gbNCname=r'\makecell{Background \\ goodness-of-fit ($\sigma$) \\ w. no corr}'
    gname=r'\makecell{Signal \\ goodness-of-fit ($\sigma$)}'
    gNCname=r'\makecell{Signal \\ goodness-of-fit ($\sigma$) \\ w. no corr}'
else:
    # Shorter names
    lname=r'\makecell{Local \\ signif. ($\sigma$)}'
    lname_twocol = r"\multicolumn{2}{c}{"+lname+r"}"
    lNCname=r'\makecell{Local \\ signif. ($\sigma$) \\ (no cor.)}'
    gbname=r'\makecell{SM \\ fit ($\sigma$)}'
    gbNCname=r'\makecell{SM \\ fit ($\sigma$) \\ (no cor)}'
    gname=r'\makecell{MSSM4 \\ fit ($\sigma$)}'
    gNCname=r'\makecell{MSSM4 \\ fit ($\sigma$) \\ (no cor)}'
new_names = [gname,gbname,lname_twocol,gNCname,gbNCname,lNCname]
for df, newname in zip(dframes,new_names):
    df.rename(   columns={use: newname}, inplace=True)

# Need to rename the degrees of freedom columns so they don't get merged between corr and no_corr datasets
DOForig = r'DOF'
DOF = r'\#SRs'
DOFNC = r'\makecell{\#SRs \\ (no cor)}'
# Renaming
for df in dframes: 
    df.rename(columns={DOForig: DOF}, inplace=True)    
# Also need to delete this column from unwanted datasets
for df in [dframes[i] for i in [1,2,4,5]]: 
    df.rename(columns={DOF: 'Null_DOF'}, inplace=True)    
# Rename no-correlation version to prevent merging
dframes[3].rename(columns={DOF: DOFNC}, inplace=True)

# Merge the dataframes
# merged = pd.merge(bg_gof,musb_0,on='experiment',how='left')
#dfs = [bg_gof,bg_gof_NC,s_gof,s_gof_NC,musb_0,musb_0_NC]
merged = reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True,how='outer'), dframes)
merged.rename(index={'Monster': 'Combined'}, inplace=True) 
# Remove underscores in experiment names, and other modifications
index = merged.index.tolist()
for name in index:
    twocol = False
    if name == "top_mass":
        newname = r"$m_{t}$"
    elif name == "alpha_s":
        newname = r"$\alpha_S$"
    elif name == "ATLAS_13TeV_MultiLEP_2Lep0Jets_36invfb":
       # Couple of special cases that we want to overlap into next column a bit
        newname = r"\multicolumn{2}{l}{ATLAS 13TeV MultiLEP 2Lep0Jets 36$\mathrm{fb}^{-1}$}"
        twocol = True
    elif name == "ATLAS_13TeV_MultiLEP_2LepPlusJets_36invfb":
        newname = r"\multicolumn{2}{l}{ATLAS 13TeV MultiLEP 2LepPlusJets 36$\mathrm{fb}^{-1}$}"
        twocol = True
    else:
        newname = name.replace("_"," ")
        newname = newname.replace("invfb","$\mathrm{fb}^{-1}$")
        newname = newname.replace("RJ3L","RJ")
    if not twocol:
        newname += " & " # Extra empty column
    merged.rename(index={name: newname}, inplace=True)

# Choose column ordering
cols = [lname_twocol,gbname,gname,DOF,lNCname,gbNCname,gNCname,DOFNC] # leaving gNCname out on purpose
merged = merged[cols]

old_names = [gNCname,gbNCname,lNCname,DOFNC]
new_names = [gname,gbname,lname,DOF]
for o,n in zip(old_names,new_names):
    merged.rename(columns={o: n}, inplace=True)

# Fix up weirdness where index column header is on a different line to other header names
#merged.columns.name = df.index.name
#merged.index.name = "Analysis"
# That doesn't work...

# Set numerical results that would render as negative or 0.0 to zero
num = merged._get_numeric_data()
num[num<0.05] = 0

print(merged)

name = "pvalue_table.tex"
print("Saving tex table as {0}".format(name))
with open(name, "w") as output_file:
    # define formatting function for columns
    def f1(x):
        if x == 0:
            s = "0"
        else:
            s = "{0:.1f}".format(x)
        return s
    def fint(x):
        return "{0:d}".format(x)
    form = [f1]*len(cols)
    form[3] = fint
    form[7] = fint
    raw_latex_string = merged.to_latex(formatters=form,escape=False)
    fixed = raw_latex_string.replace("nan","--")
    output_file.write(fixed)

# Extra rows to manually copy in
# & \multicolumn{4}{c}{ } & \multicolumn{4}{|c}{No correlations} \\ 
# \midrule


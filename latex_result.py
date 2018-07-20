"""Converting results into nice LaTeX tables"""

import pandas as pd
import pickle
from functools import reduce

pd.set_option('display.max_colwidth', -1)

with open("MSSMEW4_pvalue_results.pkl", 'rb') as pkl_file: 
    results = pickle.load(pkl_file)  

# Also get results using all signal regions, ignoring correlations when no covariance matrix is available
with open("MSSMEW4_pvalue_results_NOCORR.pkl", 'rb') as pkl_file: 
    results_NOCORR = pickle.load(pkl_file)  


# Set negative significances to zero, since that confuses people less
# Also set significances that would be rendered as 0.0 anyway to 0.
num = results.results_df._get_numeric_data()
num[num<0.05] = 0

numNC = results_NOCORR.results_df._get_numeric_data()
numNC[numNC<0.05] = 0


# Prepare latex tables for publication
# Need to rearrange them a bit, and improve labels
# Want to merge the background-only gof and musb_mu=0 results
# So pull out the underlying pandas dataframes, then add new columns to one of them
s_gof = results.query('test == "gof"').results_df
bg_gof = results.query('test == "gof_b"').results_df
musb_0 = results.query('test == "musb_mu=0"').results_df
s_gof_NC = results_NOCORR.query('test == "gof"').results_df
musb_0_NC = results_NOCORR.query('test == "musb_mu=0"').results_df

# Delete all the stuff we don't want in the result
del_list = ['test', 'a_pval', 'e_pval', 'a_sig'] # keep 'e_sig'
for df in [bg_gof,musb_0]:
    for col in del_list:
        df.pop(col)

# Do some renaming
ename='Analysis'
lname=r'\makecell{Local \\ significance ($\sigma$)}'
lNCname=r'\makecell{Local \\ significance ($\sigma$) \\ w. no corr.}'
gbname=r'\makecell{Min. global \\ significance ($\sigma$)}'
gname=r'\makecell{Signal \\ goodness-of-fit ($\sigma$)}'
gNCname=r'\makecell{Signal \\ goodness-of-fit ($\sigma$) \\ w. no corr}'
# Shorter names
#lname=r'\makecell{Local \\ signif. ($\sigma$)}'
#lNCname=r'\makecell{Local \\ signif. ($\sigma$) \\ w. no corr.}'
#gbname=r'\makecell{Min. global \\ signif. ($\sigma$)}'
#gname=r'\makecell{Signal \\ G.O.F. ($\sigma$)}'
#gNCname=r'\makecell{Signal \\ G.O.F. ($\sigma$) \\ w. no corr}'
musb_0.rename(columns={'e_sig': lname}, inplace=True)
musb_0_NC.rename(columns={'e_sig': lNCname}, inplace=True)
bg_gof.rename(columns={'e_sig': gbname}, inplace=True)
s_gof.rename(columns={'e_sig': gname}, inplace=True)
s_gof_NC.rename(columns={'e_sig': gNCname}, inplace=True)

# Merge the dataframes
# merged = pd.merge(bg_gof,musb_0,on='experiment',how='left')
dfs = [bg_gof,s_gof,musb_0,musb_0_NC,s_gof_NC]
merged = reduce(lambda left,right: pd.merge(left,right,on='experiment',how='left'), dfs)
merged.rename(columns={'experiment': ename}, inplace=True)

# Remove column with indices by setting experiment names as indices
merged.set_index(ename, inplace=True)
merged.rename(index={'Monster': 'Combined'}, inplace=True) 
# Remove underscores in experiment names, and other modifications
index = merged.index.tolist()
for name in index:
    if name == "top_mass":
        newname = r"$m_{t}$"
    elif name == "alpha_s":
        newname = r"$\alpha_S$"
    else:
        newname = name.replace("_"," ")
        newname = newname.replace("invfb","$\mathrm{fb}^{-1}$")
        newname = newname.replace("RJ3L","RJ")
    merged.rename(index={name: newname}, inplace=True)

# Choose column ordering
cols = [lname,lNCname,gbname,gname] # leaving gNCname out on purpose
merged = merged[cols]

# Fix up weirdness where index column header is on a different line to other header names
merged.columns.name = df.index.name
merged.index.name = "Analysis"

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
    raw_latex_string = merged.to_latex(formatters=[f1]*len(cols),escape=False)
    fixed = raw_latex_string.replace("nan","--")
    output_file.write(fixed)


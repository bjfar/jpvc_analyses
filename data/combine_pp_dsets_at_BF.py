"""Anders accidentally put different sets of the postprocessing results in different files.
   To make my other scripts work I need to combine them.
   
   The best fit point is located where fixedID=300000000
"""

import h5py
import numpy as np

fixedID_BF = 300000000 

file13TeV="MSSMEW_1sigma_13TeV.hdf5"
file8TeV="MSSMEW_1sigma_pp_8TeV.hdf5"
gname = "MSSMEW"

outfile="MSSMEW_LHC_BF_bothTeV.hdf5"
outf = h5py.File(outfile,'w')
outg = outf.create_group(gname)

#Copy everything at BF point from 13TeV file
hf = h5py.File(file13TeV,'r')
g = hf[gname]

# Get index of BF point
BFi = np.where(g["fixedID"][:]==fixedID_BF)[0][0]

for name, dset in g.items():
    outg.create_dataset(name, (1,), dtype='d')
    outg[name][:] = dset[BFi]

# Copy just the LHC analyses from the 8 TeV dataset (if they don't already exist in the outfile file)
hf = h5py.File(file8TeV,'r')
g = hf[gname]

# Get index of BF point
BFi = np.where(g["fixedID"][:]==fixedID_BF)[0][0]

# Copy new stuff from this point into a new HDF5 file
for name, dset in g.items():
    if name not in outg.keys(): 
        outg.create_dataset(name, (1,), dtype='d')
        outg[name][:] = dset[BFi]

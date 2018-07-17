import numpy as np

# Helper class to manage data/is_valid HDF5 dataset pairs 
class DataPair:
  def __init__(self,name,dset,dset_isvalid):
     self._name = name
     self._data  = dset
     self._valid = dset_isvalid

  def name(self):
     return self._name
 
  def data(self):
     return self._data

  def valid(self):
     return np.array(self._valid, dtype=np.bool)

  def validdata(self):
     return self.data()[self.valid()]

def get_data(in_group, hdf5_dataset_names, m=None, i=None, verbose=False):
  output = []
  for name in hdf5_dataset_names:
    if verbose: print(name,": ", len(in_group[name+"_isvalid"]))
    if m is None and i is None:
      output += [ DataPair(name, in_group[name], in_group[name+"_isvalid"]) ]
    elif m is None:
      output += [ DataPair(name, in_group[name][i], in_group[name+"_isvalid"][i]) ] 
    elif i is None:
      output += [ DataPair(name, in_group[name][m], in_group[name+"_isvalid"][m]) ]
    else:
      output += [ DataPair(name, in_group[name][m][i], in_group[name+"_isvalid"][m][i]) ]     
  return output



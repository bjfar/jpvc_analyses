"""Analysis from
   ColliderBit/src/analyses/Analysis_ATLAS_13TeV_MultiLEP_36invfb.cpp

   Likelihood construction follows ColliderBit paper 
   https://arxiv.org/pdf/1705.07919.pdf
   except that we profile instead of marginalise the nuisance parameters

   We still do the stuff regarding selection of the signal region to
   use though.
"""

import JMCtools.distributions as jtd
import JMCtools.models as jtm
import numpy as np
import scipy.stats as sps
from functools import partial
import matplotlib.pyplot as plt
from .experiment import Experiment

name = "ATLAS_13TeV_MultiLEP_36invfb" 

# Observed event counts
CMS_o = [153.,   9.,  78., 11.,  6.,  2.,  2.,  0., 11.,  4.,  3.,  9.,  0.,  
           0.,  21.,   1.,  2.,  1.,  3.,  4.]

# Background estimates (means)
CMS_b = [133.,   9.8, 68., 11.5, 2.1, 0.6, 4.1, 1.6, 4.2, 2.2, 2.8, 5.4, 1.4, 
           1.1, 21.7,  2.7, 1.6, 2.2, 1.8, 1.3]

#



"""Tools for setting up an experiment constructed from a GAMBIT
   ColliderBit analysis.
   We get information about these analyses from GAMBIT in a
   standardised format, so it is possible to mostly automate
   the construction of their pdfs.
"""

import JMCtools.distributions as jtd
import JMCtools.models as jtm
import numpy as np
import scipy as sp
import scipy.stats as sps
from functools import partial
import matplotlib.pyplot as plt
import copy
import types
from JMCtools.experiment import Experiment

# Function to map signal hypothesis into Poisson distribution parameters
# Signal+background systematics dealt with via a MULTIPLICATIVE factor
def poisson_f_mult(s, b, theta):
    l = theta*np.atleast_1d(s + b)
    #m = (l<0)
    #l[m] = 0  #Poisson cannot have negative mean
    return {'mu': l}

# Function to map signal hypothesis into Poisson distribution parameters
# Signal+background systematics dealt with via an ADDITIVE factor
def poisson_f_add(s, b, theta):
    l = np.atleast_1d(s + b + theta)
    #m = (l<0)
    #l[m] = 0  #Poisson cannot have negative mean
    #print("s:",s,", b:",b," -> l:", l)
    return {'mu': l}

# Parameter mapping function for nuisance parameter constraints
def func_nuis_corr(cov, **thetas):
    #print("in func_nuis:", thetas)
    means = np.array([thetas['theta_{0}'.format(i)] for i in range(len(thetas))])
    return {'mean': means.flatten(),
             'cov': cov}

# Gaussian constraint for additive nuisance parameter
def func_nuis_norm_add(theta,theta_std):
    #print("theta:", theta, ", theta_std:", theta_std)
    return {'loc': theta,
             'scale': theta_std}

# Log-normal constraint for multiplicative nuisance parameter
def func_nuis_lognorm_mult(theta,theta_std):
    return {'scale': theta,
            's': theta_std}

# function to join lists/arrays along last dimension
def ljoin(a,b,ma=None,mb=None):
    if ma is None: ma = np.ones(np.atleast_1d(a).shape[-1],dtype=bool) 
    if mb is None: mb = ma
    r = np.concatenate([np.atleast_1d(np.atleast_1d(a)[ma]),
                        np.atleast_1d(np.atleast_1d(b)[mb])],
                                     axis=-1)
    return r

def smooth_poisson(k,mu,loc=0):
    r =  k*np.log(mu) - mu - sp.special.gammaln(k+1)
    return r

def custpois(func,rename):
    """Construction transformed Poisson distribution,
       with logpmf function replaced by a version that
       can be evaluated for non-integer data (needed
       for Asimov calculationd"""

    mypois = jtd.TransDist(sps.poisson) # Null transformation, just to build object
    mypois.set_logpdf(smooth_poisson) # replace pdf calculation

    # Now build the transformed object
    return jtd.TransDist(mypois,func,rename)

class Analysis:
    def __init__(self,name):
        self.name = name
        # The following need to be 'manually' added to this object!
        # We use this structure because it is easier/neater to
        # autogenerate the analyses this way.
        #self.SR_names
        #self.SR_n    
        #self.SR_b    
        #self.SR_b_sys
        #self.SR_s_sys
        #self.SR_s    
        self.cov = None

    def make_experiment(self,signal=None):
        """Turn this analysis information into a jpvc experiment"""
        # Maximum likelihood estimators for 's' parameters
        # under the observed data, ignoring correlations
        # between nuisance observations. Useful for seeding
        # certain fits.
        self.s_MLE = np.array(self.SR_n) - np.array(self.SR_b)
        print('self.cov:',self.cov)
        if self.cov is None or self.cov==[]:
           e = self.make_experiment_nocov(signal)
        else:
           e = self.make_experiment_cov()
        return e        

    # Functions to provide starting guesses for parameters, tuned to each MC sample realisation
    def seeds_full_f_add(self,selected=None):
        def get_seeds_full(samples,signal):
           """Gets seeds for s and theta fits (additive nuisance)
              Gives exact MLEs for gof case!
           """
           seeds={}
           if selected is None:
               SRlist = range(self.N_SR)
               bin_samples = samples[:,0,:self.N_SR].T
               theta_samples = samples[:,0,self.N_SR:].T
           else:
               SRlist = range(self.N_SR)[selected]
               N = len(SRlist)
               bin_samples = samples[:,0,:N].T
               theta_samples = samples[:,0,N:].T
           for i,n,x in zip(SRlist,bin_samples,theta_samples):
               theta_MLE = x
               s_MLE = n - x - self.SR_b[i]
               l = s_MLE + self.SR_b[i] + theta_MLE
               # Sometimes get l < 0 predictions just due to small numerical errors. Fix these if they
               # are small enough not to matter
               mfix = (l<1e-10) & (l>-1e-10) # fix up these ones
               print('l[mfix]:', l[mfix])
               acheck = l<-1e-10 # throw error if too terrible
               if np.sum(acheck)>0:
                   raise ValueError("Negative signal prediction detected!")
               theta_MLE[mfix] = x[mfix] + 1e-10
               seeds['theta_{0}'.format(i)] = theta_MLE
               seeds['s_{0}'.format(i)] = s_MLE 
               #print('seeds for s_{0}: {1}'.format(i,s_MLE))
           return seeds
        return get_seeds_full

    def seeds_full_f_mult(self):
        def get_seeds_full(samples,full_fixed_pars):
           """Gets seeds for s and theta fits (multiplicative nuisance)"""
           seeds={}
           bin_samples = samples[:,0,:self.N_SR].T
           theta_samples = samples[:,0,self.N_SR:].T
           for i in range(self.N_SR):
              theta_MLE = theta_samples[i]
              s_MLE = bin_samples[i]/theta_MLE - self.SR_b[i]
              seeds['theta_{0}'.format(i)] = theta_MLE
              seeds['s_{0}'.format(i)] = s_MLE 
              #print('seeds for s_{0}: {1}'.format(i,s_MLE))
           return seeds
        return get_seeds_full
    
    def seeds_null_f(self,selected=None): 
        def get_seeds_null(samples,signal):
           """Gets seeds for (both additive and multiplicative) nuisance parameters fits"""
           theta_seeds={}
           if selected is None:
               SRlist = range(self.N_SR)
               theta_samples = samples[:,0,self.N_SR:].T
           else:
               SRlist = range(self.N_SR)[selected] 
               theta_samples = samples[:,0,1].T
           for i,x in zip(SRlist,theta_samples):
              theta_MLE = x
              theta_seeds['theta_{0}'.format(i)] = theta_MLE
           return theta_seeds
        return get_seeds_null

    def seeds_null_f_gof(self,selected=None,mu=1): 
        def get_seeds(samples,signal):
            """Gets seeds for additive nuisance parameters fits
               Gives exact MLEs for gof case, or when mu + s are otherwise fixed.
               Depends on signal hypothesis though, so calculation
               is a bit more involved."""
            verbose = False # For debugging
            if verbose: print("signal (seeds_null_f_gof):",signal)
            theta_seeds={}
            self.theta_both={} # For debugging, need to compare BOTH solutions to numerical MLEs.
            self.theta_dat={}
            if selected is None:
                if 2*self.N_SR!=samples.shape[-1]:
                    raise ValueError("No signal region selected, so tried to compute seeds for all signal regions; however, supplied samples do not match the shape required for all signal regions! This analysis has {0} signal regions, so we require 2*{0}={1} random variates, but we only were given {2} (samples.shape = {3})".format(self.N_SR,2*self.N_SR,samples.shape[-1],samples.shape))
                SRlist = range(self.N_SR)
                bin_samples = samples[:,0,:self.N_SR].T
                theta_samples = samples[:,0,self.N_SR:].T
            else:
                SRlist = range(self.N_SR)[selected]
                N = len(SRlist)
                bin_samples = samples[:,0,:N].T
                theta_samples = samples[:,0,N:].T
            #print("stuff:",SRlist,bin_samples,theta_samples)
            for i,n,x in zip(SRlist,bin_samples,theta_samples):
                if verbose: print("Finding MLEs for theta_{0}; total samples: {1}".format(i,len(samples)))
                if signal==None:
                   s = 0
                else:
                   s = mu*signal['s_{0}'.format(i)] #self.SR_s[i]
                b = self.SR_b[i]
                bsys = self.SR_b_sys[i]
                A = 1./bsys**2
                B = 1 + A*(s + b - x)
                C = (s+b)*(1-A*x) - n
                D = B**2 - 4*A*C
                #D<0:
                theta_MLE = np.zeros(n.shape)
                # No solutions! This is a drag, means MLE is on boundary of allowed space, not at a 'real' minima.
                # Will mess up Wilk's theorem a bit. However MLE will necessarily just be on the boundary
                # s + b + theta = 0
                # so it is easy to compute at least
                if verbose: print("s:",s)
                if verbose: print("b:",b)
                if verbose: print("D:", D.shape)
                theta_MLE[D<0] = -(s+b)
                if verbose: print("No solution count:", np.sum(D<0))
                #elif D==0:
                # One real solution:
                #print(B)
                theta_MLE[D==0] = -B[D==0]/(2*A) # Is this always positive? Not sure...
                if verbose: print("Single solution count:", np.sum(D==0))
                #elif D>0:
                # Two real solutions
                # Hmm not sure how to pick... I guess first check that they are indeed in the "allowed" region
                if verbose: print("Two solution count:", np.sum(D>0))
                r1 = (-B + np.sqrt(D))/(2*A)
                r2 = (-B - np.sqrt(D))/(2*A)
                a1 = (s+b+r1 >= 0)
                a2 = (s+b+r2 >= 0) 
                # See if just one is allowed
                ma= (D>0) & a1 & a2
                mf= (D>0) & (~a1) & (~a2)
                # m = ma | mf 
                if verbose: print("   both allowed:", np.sum(ma))
                if verbose: print("   both forbidden:", np.sum(mf))
                # #print(list(zip(n[D>0]-s-b,r1,r2)))
                # #if both allowed, pick the one closest to the MLE for theta that one would get without the nuisance constraint
                # Actually just throw error, I think this should not happen
                #if np.sum(ma)>0:
                #    for r1i, r2i in zip(r1[ma],r2[ma]):
                #        print("r1: {0}, r2: {1}, s+b: {2}".format(r1i,r2i,s+b))
                #    raise ValueError("Found multiple allowed solutions for some MLEs! I'm pretty sure this shouldn't happen so I'm calling it an error/bug in the calculation") 
                # If both are forbidden, make sure that it isn't just by some tiny amount due to numerical error.
                fix1 = ~a1 & (s+b+r1 >= -np.abs(r1)*1e-6)
                fix2 = ~a2 & (s+b+r2 >= -np.abs(r2)*1e-6)
                #m12 = (D>0) & fix1 & fix2 # Don't need to worry about this, means they are probably the same solution anyway
                m1  = (D>0) & mf & fix1
                m2  = (D>0) & mf & fix2
                #print("Corrected {0} positive solutions".format(np.sum(m1)))
                theta_MLE[m1] = -(s+b) + np.abs(r1[m1])*1e-6 # Make sure still positive after numerics
                theta_MLE[m2] = -(s+b) + np.abs(r2[m2])*1e-6
                if np.sum(mf & ~(fix1 | fix2) > 0):
                    for r1i, r2i in zip(r1[mf],r2[mf]):
                        print("r1: {0}, r2: {1}, s+b: {2}".format(r1i,r2i,s+b))
                    raise ValueError("Found both solutions forbidden (and unfixable) for some MLEs! I'm pretty sure this shouldn't happen so I'm calling it an error/bug in the calculation") 
                #d1 = (r1 - (n - s - b))**2
                #d2 = (r2 - (n - s - b))**2
                # Or the closest to the MLE using a normal approximation for the Poisson? And fixing variance with theta=0...
                MLE_norm = (x*A + (n-s-b)/(s+b)) / (A + 1./(s+b))
                d1 = (r1 - MLE_norm)**2
                d2 = (r2 - MLE_norm)**2
                # Use this if both solutions are allowed.
                m1 = (D>0) & a1 & a2 & (d1<d2)
                m2 = (D>0) & a1 & a2 & (d1>=d2)
                if verbose: print("   both roots allowed, but positive is closer to gaussian estimate:", np.sum(m1))
                if verbose: print("   both roots allowed, but negative is closer to gaussian estimate:", np.sum(m2))
                theta_MLE[m1] = r1[m1]
                theta_MLE[m2] = r1[m2]
                # or just smallest in magnitude?
                #d1 = r1**2
                #d2 = r2**2
                # The most proper thing is just to see which solution gives the best likelihood (or minimum log-likelihood)
                #d1 = -sps.poisson.logpmf(n,s+b+r1) - sps.norm.logpdf(x,r1)
                #d2 = -sps.poisson.logpmf(n,s+b+r2) - sps.norm.logpdf(x,r2)
                #d1[~np.isfinite(d1)] = 1e99 # bad likelihood
                #d2[~np.isfinite(d2)] = 1e99
                #theta_MLE[(D>0) & (d1<d2)] = r1[(D>0) & (d1<d2)]
                #theta_MLE[(D>0) & (d1>=d2)] = r2[(D>0) & (d1>=d2)]
                #print("r1", "r2", "MLE_guess")
                #for j in range(np.sum(D>0)):
                #    print(r1[j], r2[j], MLE_norm[j])
                #theta_MLE[D>0] = MLE_norm[D>0] # test...
                # Seems like the positive root is always the right one for some reason...
                #theta_MLE[D>0] = r1[D>0] 
                # # If both solutions are forbidden then something weird is probably going on, but we'll pick one anyway
                # # and let the parameter mapping sort it out.
                m = (D>0) & (a1) & (~a2)
                if verbose: print("   only positive root allowed:", np.sum(m))
                theta_MLE[m] = r1[m]
                m = (D>0) & (~a1) & (a2)
                if verbose: print("   only negative root allowed:", np.sum(m))
                theta_MLE[m] = r2[m]
                # Save some extra info for diagnositic functions
                self.theta_both['theta_{0}'.format(i)] = [copy.copy(theta_MLE), copy.copy(theta_MLE), MLE_norm, s + b] # Make sure to copy...
                self.theta_both['theta_{0}'.format(i)][0][D>0] = r1[D>0]
                self.theta_both['theta_{0}'.format(i)][1][D>0] = r2[D>0]
                self.theta_dat['theta_{0}'.format(i)] = n,x,s,b,bsys
                # Output!
                if verbose: print("MLE theta_{0}:".format(i), theta_MLE) 
                # Sanity check; are all these seeds part of the allowed parameter space?
                l = s + b + theta_MLE
                afix = (l<=1e-10) & (l>-1e-10) # really small values can cause problems too.
                theta_MLE[afix] = theta_MLE[afix] + 1e-10
                acheck = l<-1e-10 # Raise error if too far negative to correct without problems
                if np.sum(acheck)>0:
                     for soli, r1i, r2i in zip(theta_MLE[acheck],r1[acheck],r2[acheck]):
                        print("MLE: {0}, r1: {1}, r2: {2}, s+b: {3}, ? {4}".format(soli,r1i,r2i,s+b,s+b+r1i >= -np.abs(r1i)*1e-6))
                     raise ValueError("Computed {0} forbidden seeds (from {1} samples)! There is therefore a bug in the seed calculations".format(np.sum(acheck),len(samples)))
                theta_seeds['theta_{0}'.format(i)] = theta_MLE
            return theta_seeds
        return get_seeds

    def seeds_null_f_cov(self): 
        def get_seeds_null_cov(samples,signal):
            """Gets seeds for additive nuisance parameters fits with
               covariance matrix. I hope they are exact, need to test...
               Nope, this is just incorrect."""
            theta_seeds={}
            bin_samples = samples[:,0,:self.N_SR].T
            theta_samples = samples[:,0,self.N_SR:].T
            # Need inverse of covariance matrix
            V = np.linalg.inv(self.cov)
            for k in range(self.N_SR):
                s = signal['s_{0}'.format(k)] 
                b = self.SR_b[k]
                Uk = V[k][k] * (s+b)
                A = 1
                B = 4*Uk + 1
                Xi = np.sum(V[k,:,np.newaxis] * theta_samples,axis=0)
                C = 8*Uk**2 - 4*V[k][k]*bin_samples[k] - 4*Uk*Xi 
                Jp = (-B + np.sqrt(B**2 - 4*A*C)) / 2.
                Jm = (-B - np.sqrt(B**2 - 4*A*C)) / 2.
                print('Uk:', Uk)
                print('Xi:', Xi)
                print('C:', C)
                print('B**2 - 4*A*C:', B**2 - 4*A*C)
                print('Jp:', Jp)
                print('Jm:', Jm)           
                theta_seeds['theta_{0}'.format(k)] = (Uk + Jp/2.) / V[k][k] 
            print(theta_seeds)
            return theta_seeds
        return get_seeds_null_cov

 
    def make_get_asimov_nocov(self,selected=None):
        """Returns a function which returns the Asimov data for this analysis
           (with respect to signal-strength 'mu' parameters)
           Works for both correlated and uncorrelated cases (so long as nuisance
           parameter true values are zero in the simulated data).
        """
        def get_asimov_data(mu,signal=None):
            if mu==0: # background-only case
                nA = ljoin(self.SR_b,np.zeros(self.N_SR),selected)
            else: # signal + background case; usually only want mu=1 but I allowed the more general case
                s = [signal['s_{0}'.format(i)] for i in range(self.N_SR)]
                nA = ljoin(np.array(self.SR_b) + mu*np.array(s),np.zeros(self.N_SR),selected)
            # We also need the MLEs of nuisance parameters under this Asimov data
            nuis_MLEs = {'theta_{0}'.format(i): 0 for i in range(self.N_SR)}
            return nA, nuis_MLEs
        return get_asimov_data

    def make_experiment_cov(self):
        # Create the transformed pdf functions
        # Also requires some parameter renaming since we use the
        # same underlying function repeatedly
        poisson_part = [custpois(partial(poisson_f_add,b=self.SR_b[i]),
                               ['s_{0} -> s'.format(i), 
                                'theta_{0} -> theta'.format(i)])
                         for i in range(self.N_SR)]
        corr_dist = jtd.TransDist(sps.multivariate_normal,partial(func_nuis_corr,cov=self.cov),
                       func_args=["theta_{0}".format(i) for i in range(self.N_SR)])
        correlations = [(corr_dist,self.N_SR)]

        # Create the joint PDF object
        joint = jtd.JointDist(poisson_part + correlations)
         
        # Set options for parameter fitting
        theta_opt  = {'theta_{0}'.format(i) : 0 for i in range(self.N_SR)}
        theta_opt2 = {'error_theta_{0}'.format(i) : 0.1*np.sqrt(self.cov[i][i]) for i in range(self.N_SR)} # Get good step sizes from covariance matrix
        s_opt  = {'s_{0}'.format(i): 0 for i in range(self.N_SR)} # Maybe zero is a good starting guess? Should use seeds that guess based on data.
        s_opt2 = {'error_s_{0}'.format(i) :  0.1*np.sqrt(self.cov[i][i]) for i in range(self.N_SR)} # Get good step sizes from covariance matrix.
        s_options = {**s_opt, **s_opt2}
        
        nuis_options = {**theta_opt, **theta_opt2}
        general_options = {**s_options, **nuis_options}

        # Full observed data list, included observed values of nuisance measurements
        observed_data = ljoin(self.SR_n, np.zeros(self.N_SR))

        # Define the experiment object and options for fitting during statistical tests
        e = Experiment(self.name,joint,observed_data,DOF=self.N_SR)
    
        e.define_gof_test(null_options=nuis_options,
                          full_options=general_options,
                          null_seeds=(self.seeds_null_f_gof(), False), # Seeds NOT exact with covariance matrix! Just testing.
                          full_seeds=(self.seeds_full_f_add(), False), 
                          diagnostics=[self.make_dfull(s_opt,theta_opt),
                                       self.make_dnull(theta_opt),
                          ])
        
        e.define_mu_test(null_options=nuis_options,
                         null_seeds=(self.seeds_null_f_gof(), False),
                         scale_with_mu=list(s_opt.keys()),
                         )

        e.define_musb_test(null_options=nuis_options,
                           mu1_seeds=(self.seeds_null_f_gof(mu=1), False), # naming a bit odd, but these are the mu=1 seeds
                           mu0_seeds=(self.seeds_null_f_gof(mu=0), False), # " "   mu=0
                           scale_with_mu=list(s_opt.keys()),
                           asimov=self.make_get_asimov_nocov() # pretty sure Asimov data is the same regardless of correlations.
                           )
        return e


    def make_experiment_nocov(self,signal):
        # Create the transformed pdf functions
        # Also requires some parameter renaming since we use the
        # same underlying function repeatedly
        # poisson_part_mult = [jtd.TransDist(sps.poisson,partial(poisson_f_mult,b=self.SR_b[i]),
        #                        ['s_{0} -> s'.format(i), 
        #                         'theta_{0} -> theta'.format(i)])
        #                  for i in range(self.N_SR)]

        poisson_part_add = [custpois(partial(poisson_f_add,b=self.SR_b[i]),
                               ['s_{0} -> s'.format(i), 
                                'theta_{0} -> theta'.format(i)])
                         for i in range(self.N_SR)]

        # Using lognormal constraint on multiplicative systematic parameter
        # sys_dist_mult = [jtd.TransDist(sps.lognorm,
        #                           partial(func_nuis_lognorm_mult,
        #                                   theta_std=self.SR_b_sys[i]/self.SR_b[i]),
        #                           ['theta_{0} -> theta'.format(i)])
        #               for i in range(self.N_SR)]

        # Using normal constaint on additive systematic parameter
        sys_dist_add = [jtd.TransDist(sps.norm,
                                  partial(func_nuis_norm_add,
                                          theta_std=self.SR_b_sys[i]),
                                  ['theta_{0} -> theta'.format(i)])
                      for i in range(self.N_SR)]

        # Median data under background-only hypothesis
        expected_data = ljoin(np.round(self.SR_b), np.zeros(self.N_SR))
        expected_data = expected_data[np.newaxis,np.newaxis,:] # Add required extra axes.

        #print("fractional systematic uncertainties:")
        #print([self.SR_b_sys[i]/self.SR_b[i] for i in range(self.N_SR)])
        #quit()

        # This next part is a little tricky. We DON'T know the correlations
        # between signal regions here, so we follow the method used in
        # ColliderBit and choose just one signal region to use in our test,
        # by picking, in advance, the region with the best sensitivity to
        # the signal that we are interested in.
        # That is, the signal region with the highest value of
        # Delta LogL = LogL(n=b|s,b) - LogL(n=b|s=0,b)
        # is selected.
        #
        # So, we need to compute this for all signal regions.
        seedf = self.seeds_null_f_gof()
        seedb = seedf(expected_data,signal) # null hypothesis fits depend on signal parameters
        zero_signal = {'s_{0}'.format(i): 0 for i in range(self.N_SR)}
        seed  = seedf(expected_data,zero_signal)
        LLR = []
        for i in range(self.N_SR):
            model = jtm.ParameterModel([poisson_part_add[i]]+[sys_dist_add[i]])
  
            odatai = np.array([np.round(self.SR_b[i])]+[0]) # median expected background-only data
            si = 's_{0}'.format(i)
            ti = 'theta_{0}'.format(i)
            parsb = {ti: seedb[ti], **zero_signal}
            pars  = {ti: seed[ti], **signal}

            Lmaxb = model.logpdf(parsb,odatai)
            Lmax  = model.logpdf(pars,odatai)

            LLR += [-2 * (Lmax - Lmaxb)]
           
        # Select region with largest expected (background-only) LLR for this signal
        # (Note, if input signal is in fact zero, LLR will be zero for all signal regions, and
        # signal region zero will always get chosen)
        selected = slice(np.argmax(LLR),np.argmax(LLR)+1) # keep slice format for generality

        # Just for curiosities sake, let's disable the signal region selection and treat them all as independent:
        #selected = slice(0,self.N_SR)

        print("Selected signal region {0} ({1}) in analysis {2}".format(selected,self.SR_names[selected],self.name))
        submodels = poisson_part_add[selected] + sys_dist_add[selected]

        # Create the joint PDF object
        #joint = jtd.JointDist(poisson_part_mult + sys_dist_mult)
        joint = jtd.JointDist(submodels) 

        sel_i = range(self.N_SR)[selected]
        theta_opt  = {'theta_{0}'.format(i) : 0 for i in sel_i} # additive
        theta_opt2 = {'error_theta_{0}'.format(i) : 1.*self.SR_b_sys[i] for i in sel_i} # Get good step sizes from systematic error estimate
        s_opt  = {'s_{0}'.format(i): 0 for i in sel_i} # Maybe zero is a good starting guess? Should use seeds that guess based on data.
        s_opt2 = {'error_s_{0}'.format(i) :  0.1*self.SR_b_sys[i] for i in sel_i} # Get good step sizes from systematic error estimate
        s_options = {**s_opt, **s_opt2}
       
        nuis_options = {**theta_opt, **theta_opt2} #, 'print_level':1}
        general_options = {**s_options, **nuis_options}  

        #print("nuis_options   :", nuis_options)
        #print("general_options:", general_options)

        # # Set options for parameter fitting
        # #theta_opt  = {'theta_{0}'.format(i) : 1 for i in range(self.N_SR)} # multiplicative
        # theta_opt  = {'theta_{0}'.format(i) : 0 for i in range(self.N_SR)} # additive
        # theta_opt2 = {'error_theta_{0}'.format(i) : 1.*self.SR_b_sys[i] for i in range(self.N_SR)} # Get good step sizes from systematic error estimate
        # s_opt  = {'s_{0}'.format(i): 0 for i in range(self.N_SR)} # Maybe zero is a good starting guess? Should use seeds that guess based on data.
        # s_opt2 = {'error_s_{0}'.format(i) :  0.1*self.SR_b_sys[i] for i in range(self.N_SR)} # Get good step sizes from systematic error estimate
        # s_options = {**s_opt, **s_opt2}
       
        # nuis_options = {**theta_opt, **theta_opt2} #, 'print_level':1}
        # general_options = {**s_options, **nuis_options}

        # print("Setup for experiment {0}".format(self.name))
        # #print("general_options:", general_options)
        # #print("s_MLE:", self.s_MLE)
        # #print("N_SR:", self.N_SR)
        # #print("observed_data:", observed_data.shape)
        # oseed = self.seeds_full_f_mult()(np.array(observed_data)[np.newaxis,np.newaxis,:])
        # print("parameter, MLE, data, seed")
        # for i in range(self.N_SR):
        #     par = "s_{0}".format(i)
        #     print("{0}, {1}, {2}, {3}".format(par, self.s_MLE[i], observed_data[i], oseed[par]))
        # for i in range(self.N_SR):
        #     par = "theta_{0}".format(i)
        #     print("{0}, {1}, {2}, {3}".format(par, 1, observed_data[i+self.N_SR], oseed[par]))
        # quit()

        # Define the experiment object and options for fitting during statistical tests
        print(selected)
        print(np.array(self.SR_n)[selected])
        print(np.zeros(self.N_SR)[selected])
        odata = ljoin(np.round(self.SR_n), np.zeros(self.N_SR), selected) 
        e = Experiment(self.name,joint,odata,DOF=len(sel_i))
         
        e.define_gof_test(null_options=nuis_options,
                          full_options=general_options,
                          null_seeds=(self.seeds_null_f_gof(selected), True),
                          full_seeds=(self.seeds_full_f_add(selected), True), # Extra flag indicates that the "seeds" are actually the analytically exact MLEs, so no numerical minimisation needed
                          diagnostics=[self.make_dfull(s_opt,theta_opt,selected),
                                       self.make_dnull(theta_opt,selected),
                          ])
                          #             self.make_seedcheck(),
                          #             self.make_checkpdf()]
                          #)
        
        e.define_mu_test(null_options=nuis_options,
                         null_seeds=self.seeds_null_f_gof(selected),
                         scale_with_mu=['s_{0}'.format(i) for i in sel_i],
                         )

        e.define_musb_test(null_options=nuis_options,
                           mu1_seeds=(self.seeds_null_f_gof(selected,mu=1), True), # naming a bit odd, but these are the mu=1 seeds
                           mu0_seeds=(self.seeds_null_f_gof(selected,mu=0), True), # " "   mu=0
                           scale_with_mu=['s_{0}'.format(i) for i in sel_i],
                           asimov=self.make_get_asimov_nocov(selected)
                           )

        # Just check that pdf calculation gives expected answer:
        # pars = {**s_opt,**theta_opt}
        # x = np.zeros(self.N_SR)
        # logpdf = e.general_model.logpdf(pars,e.observed_data)
        # expected_logpdf = [sps.poisson.logpmf(self.SR_n[i],self.SR_b[i]+pars['s_{0}'.format(i)]+pars['theta_{0}'.format(i)]) for i in range(self.N_SR)] \
        #                   + [sps.norm.logpdf(x[i],loc=pars['theta_{0}'.format(i)],scale=self.SR_b_sys[i]) for i in range(self.N_SR)]
        # print('logpdf         :',logpdf)
        # print('expected logpdf:', np.sum(expected_logpdf))

        # print("Components:")
        # for l, el in zip(e.general_model.logpdf_list(pars,e.observed_data), expected_logpdf):
        #     print('   logpdf:{0},  exp:{1}'.format(l[0][0],el))

        return e

    def make_dfull(self,s_opt,theta_opt,selected=None):
        # Can define extra calculations to be done or plots to be created using the fit
        # results, to help diagnose any problems with the fits. 
        def dfull(e, Lmax0, pmax0, seeds0, Lmax, pmax, seeds, samples):
            # Plot distribution of fit values against their
            # true values under the null hypothesis. Make sure
            # this makes sense.
         
            expected = {**s_opt,**theta_opt}
        
            fig = plt.figure(figsize=(2*self.N_SR,6))
            N = len(pmax.keys())
           
            if selected is None:
               SRlist = range(self.N_SR)
            else:
               SRlist = range(self.N_SR)[selected]
            NSR = len(SRlist)
            for i in range(2*NSR):
                if i % 2 == 0:
                   pos = i//2 + 1
                   key = 's_{0}'.format(SRlist[i//2])
                elif i % 2 == 1:
                   pos = i//2 + 1 + NSR
                   key = 'theta_{0}'.format(SRlist[i//2])
                val = pmax[key] 
                val = val[np.isfinite(val)] # remove nans from failed fits
                #val = val[val>=0] # remove non-positive MLEs, these can't be log'd
                n, bins = np.histogram(val, normed=True)
                ax = fig.add_subplot(2,NSR,pos)
                ax.plot(bins[:-1],n,drawstyle='steps-post',label="")
                ax.set_title(key)
                trueval = expected[key]
                ax.axvline(trueval,lw=2,c='k')
            plt.tight_layout()
            fig.savefig('{0}_diagnostic_full.png'.format(e.name))
            plt.close(fig)

        return dfull            

    def make_dnull(self,theta_opt,selected=None):
        def dnull(e, Lmax0, pmax0, seeds0, Lmax, pmax, seeds, samples):
            # Plot distribution of fit values against their
            # true values under the null hypothesis. Make sure
            # this makes sense.
         
            expected = {**theta_opt}
        
            fig = plt.figure(figsize=(2*self.N_SR,3))
            N = len(pmax.keys())

            if selected is None:
               SRlist = range(self.N_SR)
            else:
               SRlist = range(self.N_SR)[selected]
            NSR = len(SRlist)
            for i,SR in enumerate(SRlist):
                key = 'theta_{0}'.format(SR)
                pos = i+1
                val = np.atleast_1d(pmax0[key])
                val = val[np.isfinite(val)] # remove nans from failed fits
                #val = val[val>=0] # remove non-positive MLEs, these can't be log'd 
                #print(key, val)
                n, bins = np.histogram(val, normed=True)
                ax = fig.add_subplot(1,NSR,pos)
                ax.plot(bins[:-1],n,drawstyle='steps-post',label="")
                ax.set_title(key)
                trueval = expected[key]
                ax.axvline(trueval,lw=2,c='k')
            plt.tight_layout()
            fig.savefig('{0}_diagnostic_null.png'.format(e.name))
            plt.close(fig)
        return dnull

    def make_seedcheck(self):
        def dseed(e, Lmax0, pmax0, seeds0, Lmax, pmax, seeds, samples):
            # Compare seeds to MLEs
            print("Comparing seeds to final MLEs for null hypothesis")
            solnames = ['seed','MLE ','r1  ','r2  ','rn  ']
            relogpdf = {n:np.zeros(len(samples)) for n in solnames}
            for key in seeds0.keys():
                print("  {0}:".format(key))
                #print(seeds0[key])
                #print(pmax0[key])
                for i,(seed, MLE) in enumerate(zip(np.atleast_1d(seeds0[key]), np.atleast_1d(pmax0[key]))):
                    r1, r2, rn, sb = self.theta_both[key]
                    if seed!=MLE:
                        #print("     seed:{0},  MLE:{1}".format(seed,MLE))
                        print("     seed:{0},  MLE:{1}, (r1={2}, r2={3}, rn={4}, s+b={5})".format(seed,MLE,r1[i],r2[i],rn[i],sb))
                    # Check likelihood values for each of these, make sure MLE is really the best fit!
                    n,x,s,b,b_sys = self.theta_dat['theta_{0}'.format(i)] 
                    for k,val in enumerate([seed,MLE,r1[i],r2[i],rn[i]]):
                        logpdf = sps.poisson.logpmf(n[i],mu=s+b+val) + sps.norm.logpdf(x[i],loc=val,scale=b_sys)
                        relogpdf[solnames[k]][i] += logpdf
                        print("       -2*logpdf: {0}, s+b+val: {1}, (running tot:{2})".format(-2*logpdf,s+b+val, -2*relogpdf[solnames[k]][i]))
 
            print("Combined likelihood using various parameters:")
            for i in range(len(samples)):
                print("Sample {0}".format(i))
                pmaxi = {par:val[i] for par,val in pmax0.items()}
                pars = {name:{**pmaxi} for name in solnames}
                for key in seeds0.keys(): 
                    seed = seeds0[key][i]
                    MLE  = pmax0[key][i]
                    r1, r2, rn, sb = self.theta_both[key]
                    pars['seed'][key] = seed
                    pars['MLE '][key] = MLE
                    pars['r1  '][key] = r1[i]
                    pars['r2  '][key] = r2[i]
                    pars['rn  '][key] = rn[i]
                minl = 1e99
                minn = ""
                for name in pars.keys(): 
                    logpdf = e.general_model.logpdf(pars[name],x=samples[i])
                    if -2*logpdf < minl: 
                        minn = name
                        minl = -2*logpdf
                    print('   -2*logpdf with {0}: {1} -- recalc: {2}'.format(name,-2*logpdf,-2*relogpdf[name][i]))
                print('   Minima was at {0} value'.format(minn))
            print()
            print("Comparing seeds to final MLEs for full hypothesis")
            for key in seeds.keys():
                print("  {0}:".format(key))
                for seed, MLE in zip(np.atleast_1d(seeds[key]), np.atleast_1d(pmax[key])):
                    if seed!=MLE:
                        print("     seed:{0},  MLE:{1}".format(seed,MLE))
        return dseed

    def make_checkpdf(self):
        def dpdf(e, Lmax0, pmax0, seeds0, Lmax, pmax, seeds, samples):
            # Checking the pdf calculation at the best fit points
            print("Manually recalculating logpdf values")
            solnames = ['seed','MLE ','r1  ','r2  ','rn  ']
            for i in range(len(samples)):
                print("Sample {0}".format(i))
                pmaxi = {par:val[i] for par,val in pmax0.items()}
                pars = {name:{**pmaxi} for name in solnames}
                for key in seeds0.keys(): 
                    seed = seeds0[key][i]
                    MLE  = pmax0[key][i]
                    r1, r2, rn, sb = self.theta_both[key]
                    pars['seed'][key] = seed
                    pars['MLE '][key] = MLE
                    pars['r1  '][key] = r1[i]
                    pars['r2  '][key] = r2[i]
                    pars['rn  '][key] = rn[i]
                for name in pars.keys():
                    print("Using {0} as parameter value".format(name)) 
                    logpdf = e.general_model.logpdf_list(pars[name],x=samples[i])
                    # recalculate pdf components manually
                    relogpdf = []
                    for j in range(self.N_SR):
                        t = 'theta_{0}'.format(j)
                        val = pars[name][t]
                        n,x,s,b,b_sys = self.theta_dat[t]
                        relogpdf += [sps.poisson.logpmf(n[i],mu=s+b+val)]
                    for j in range(self.N_SR):
                        t = 'theta_{0}'.format(j)
                        val = pars[name][t]
                        n,x,s,b,b_sys = self.theta_dat[t]
                        relogpdf += [sps.norm.logpdf(x[i],loc=val,scale=b_sys)]
                    for k, (c, ce) in enumerate(zip(logpdf, relogpdf)):
                        print('   component {0}:  pdf: {1},  re-pdf: {2}'.format(k,c,ce))
        return dpdf


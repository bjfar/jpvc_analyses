"""Script to recompute the likelihood for the CMS analyses that
   use correlation information (was a bug in the GAMBIT implementation)"""

import numpy as np
import scipy.stats as sps
import time
from copy import copy
import h5py
import pickle
import datatools as dt
import experiments.CBit_LHC_python as LHC
import mpi4py.MPI as MPI

def compute_marg_logl(n,l,cov,M,tol=0.01,rank=0):
    N_SR = len(n)
    t0 = time.time()
    converged = False
    like = 0
    like_prev = 0
    like_sum = 0
    like_sum_prev = 0
    Mnext = M
    Mtot = 0
    first = True
    while not converged:
        X = sps.multivariate_normal.rvs(mean = np.zeros(N_SR), cov=cov, size=Mnext)
        Mtot += Mnext
    
        # Do averaging
        logl = 0
        lm = np.atleast_1d(l)[:,np.newaxis] + X.T
        # Dump batch of random numbers to file for checking
        #np.savetxt("py_l.txt",lm.T)
        #quit()

        lm[lm<=1e-3] = 1e-3
        # Do what ColliderBit does, and set tiny rates to 1e-3
        newlogl = np.sum(np.atleast_1d(n)[:,np.newaxis]*np.log(lm) - lm,axis=0) # sum over signal regions
        # Dump result for checking
        np.savetxt("py_l.txt",newlogl)
        quit()

        logl += newlogl
        # Make sure to do average over likelihood, not log-likelihood...
        likes  = np.exp(logl)
        #newlike = np.mean(likes,axis=-1)
        # Combine with previous estimate
        #like = (like_prev + newlike) / 2.
        like_sum  = np.sum(likes)
        like      = like_sum / Mnext

        if first:
            first = False
            like_sum_prev = like_sum
            pass # need two batches before we start doubling
        else:
            Mnext *= 2 # Double size of next batch (ensures that we can combine batchs with equal weights)
            if Mnext > 1e7:
                converged = True # Stop if sample size getting too big
    
        # Check convergence
        pdiff = np.abs((like - like_prev) / like)
        #print("pdiff:", pdiff, "Mtot:", Mtot)
        if pdiff < tol:
            converged = True
        like_sum_prev += like_sum
        like_prev = like_sum_prev / Mtot
                
    loglike = np.log(like_prev)
    print("(rank {0}) time taken: {1} seconds (total samples: {2}, accuracy: {3}, result: {4})".format(rank, time.time() - t0, Mtot, pdiff, loglike))
    return loglike 

if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    dfile = "data/MSSMEW_pp_final.hdf5"
    
    gname = "MSSMEW"
    f = h5py.File(dfile,'r')
    g = f[gname]
    
    # Analyses whose likelihoods we need to recompute
    recompute_for = [
     "CMS_13TeV_2LEPsoft_36invfb"
    ,"CMS_13TeV_2OSLEP_36invfb"
    ]
    
    # Names of datasets containing signal predictions and their MC uncertainties
    signals = {}
    s_uncert = {}
    loglikes = {}
    new_loglikes = {}
    analyses = [a for a in LHC.analyses if a.name in recompute_for]
    for a in analyses:
        signals[a.name]  = ["#LHC_signals @ColliderBit::calc_LHC_signals::{0}__{1}__signal".format(a.name,SR) for SR in a.SR_names]
        s_uncert[a.name] = ["#LHC_signals @ColliderBit::calc_LHC_signals::{0}__{1}__signal_uncert".format(a.name,SR) for SR in a.SR_names]
        loglikes[a.name] = "#LHC_LogLike_per_SR @ColliderBit::get_LHC_LogLike_per_SR::{0}__combined_LogLike".format(a.name)
    
    # Combined likelihoods (will need to re-weight these)
    LHC_loglike = "#LHC_Combined_LogLike @ColliderBit::calc_combined_LHC_LogLike"
    LogLike = "LogLike"
    
    llike = dt.get_data(g, [LogLike])[0] #, m, i)
 
    # maximum likelihood point
    maxi = np.nanargmax(llike.data())

    N = len(llike.data())
    #N = 1 # testing
    
    # Get MSSM4 parameters so that we can cross-check ColliderBit calculations for same point
    par_dsets=[
    "#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::Ad_3"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::Ae_3"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::Au_3"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::M1"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::M2"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::M3"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::Qin"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::TanBeta"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::mA"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::ml2"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::mq2"
    ,"#MSSM11atQ_mA_parameters @MSSM11atQ_mA::primary_parameters::mu"
    ]
    Ad3,Ae3,Au3,M1,M2,M3,Qin,TanBeta,mA,ml2,mq2,mu = dt.get_data(g, par_dsets)

    def printpoint(i):
        print("Ad3:",Ad3.data()[i])
        print("Ae3:",Ae3.data()[i])
        print("Au3:",Au3.data()[i])
        print("M3 :",M3 .data()[i])
        print("Qin:",Qin.data()[i])
        print("mA :",mA .data()[i])
        print("ml2:",ml2.data()[i])
        print("mq2:",mq2.data()[i])
        print("M1 :",M1 .data()[i])
        print("M2 :",M2 .data()[i])
        print("mu :",mu .data()[i])
        print("TanBeta:",TanBeta.data()[i])

    print("Number of parameter points to process:",N)
    
    # Number of MC samples to use for marginalisation
    M = int(1e5) # Starting number of samples to use for marginalisation
    
    # Actually I guess we should just compute this when we do the final insertion back into the HDF5 file
    # I.e. separate to the processing job.
    # First we need the marginal SM likelihoods, i.e. with signal=0
    # logl_b = {}
    # X = {}
    # for a in analyses:
    #     logl_b[a.name] = compute_marg_logl(a.SR_n,a.SR_b,a.cov,M)
    # 
    # print("logl_b:",logl_b)
    
    # Now for the MSSM4 likelihoods
     
    # Results storage
    if rank==0:
        final_loglikes = np.zeros((N,len(analyses)))
    
    batch_size = 1 #50
    
    # Need to do the MPI in a smart way.
    # I think best to have rank 0 manage file IO,
    # and distribute work to worker nodes on request
    # Rank 0 will check for messages after every
    # likelihood it processes, whilst other
    # ranks only communicate back when they need
    # more work to do.
    
    globali = maxi #HACK 0 # index of next point to be processed
    pos = [None for x in range(size)]
    process_finished = [False for x in range(size)]
    
    def next_batch(i):
        # DEBUGGING ONLY
        global globali
        printpoint(globali) # Only use when mpisize=1 and batch_size=1
        # Get next batch of work for process i
        start = globali
        if start > N:
            this_batch_size = 0
        else:
            if i==0:
                new_batch_size = 1 # Rank 0 needs to check for work requests often
            else:
                new_batch_size = batch_size
            end = globali + new_batch_size
            if end > N:
                end = N
            this_batch_size = end - start
            chunk = slice(start,end) 
            pos[i] = chunk # Save location data for this slice
            globali = end
    
            # Read in data for this chunk
            s = []
            serr = []
            for k,a in enumerate(analyses):
                #print("A: time elapsed: {0} seconds".format(time.time() - starttime))
                sdata   = dt.get_data(g, signals[a.name], i=chunk)
                suncert = dt.get_data(g, s_uncert[a.name], i=chunk)
                s    += [np.array([x.data() for x in sdata])]
                serr += [np.array([u.data() for u in suncert])]
        if this_batch_size==0: 
            this_batch_size = -1 # Send done signal
            s = None
            serr = None
            process_finished[i] = True 
        return this_batch_size, s, serr
    
    def save_done_batch(batch,i):
        if pos[i] is not None:
            # Insert batch of work done by process i into results array
            start = pos[i].start
            end = pos[i].stop
            length = end-start
            if pos[i] is None:
                raise ValueError("No slice data found for process {0}! It must not have been saved correctly!".format(i))
            final_loglikes[start:end] = batch[:length]
        else:
            # We didn't give them any work to do, so this is just a request for their first lot of work
            pass
    
    def done_length():
        # Find the highest index below which everything is finished computing
        i = 0
        for slicepos in pos:
            if i==0 and slicepos is not None:
                i = slicepos.start
            elif slicepos is not None and slicepos.start < i:
                i = slicepos.start
        return i
    
    job_finished = False
    my_batch = None
    my_batch_size = 0
    new_loglikes = np.zeros((batch_size,len(analyses)))
    starttime = time.time()
    t0 = -1e99
    while not job_finished:
        if my_batch is not None and my_batch[0] != -1:
            my_batch_size, s_batch, s_err_batch = my_batch
            for i in range(my_batch_size):
                # Cycle through analyses
                for k,a in enumerate(analyses):
                    s = s_batch[k][:,i]
                    serr = s_err_batch[k][:,i]
                    cov = np.array(a.cov) + np.diag(serr**2)
                    l = s + np.array(a.SR_b)
                    new_loglikes[i,k] = compute_marg_logl(a.SR_n,l,cov,M)
    
        if rank==0:
            # Save last batch of results
            newdata = False
            if my_batch is not None and my_batch[0]!=-1:
                save_done_batch(new_loglikes,0)
                newdata = True
           
            # Check for work requests from other processes
            status = MPI.Status()
            print("Master is checking if other processes need more work...")
            while comm.Iprobe(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status):
                rnk = status.Get_source()
                # Keep looping until all work requests filled
                data = comm.recv(source=rnk)
                newdata = True
                print("Master received new results from rank {0}".format(rnk))
                save_done_batch(data,rnk)
                # Send new batch back to rnk
                newbatch = next_batch(rnk)
                comm.isend(newbatch, dest=rnk)
                if newbatch[0] == -1:
                    print("Master has informed rank {0} that there are no more points to process".format(rnk))
                else:
                    print("Master sent a new batch of {0} points to rank {1}".format(newbatch[0],rnk))
    
            # Let's save some interim results, why not
            if newdata:
                print("Pickling results into new_CMS_likes.pkl")
                with open("new_CMS_likes_partial.pkl", 'wb') as pkl_file: 
                    pickle.dump(final_loglikes[:done_length()],pkl_file)
     
            if time.time() - t0 < 2:
                # Wait two seconds
                print("Master is waiting 2 seconds...") 
                time.sleep(2)
            t0 = time.time()
    
            # Get data for self
            my_batch = next_batch(0)
    
            if np.all(process_finished):
                job_finished = True
    
        else:
            # Finished batch; send it back to the master and ask for more work
            print("Rank {0} is sending {1} results back to master".format(rank,my_batch_size))
            comm.send(new_loglikes, dest=0)
            print("Rank {0} is waiting to receive more work from master...".format(rank))
            my_batch = comm.recv(source=0)
            if my_batch[0] == -1:
                job_finished = True
                print("Rank {0} has received 'job done' signal from master and will prepare for MPI_finalise".format(rank))
            else:
                print("Rank {0} has received {1} new points from master".format(rank,my_batch[0]))
    
    # Pickle results
    print("Rank {0} finished!".format(rank))
    endtime = time.time()
    print("\nRank {0} took {1:.0f} seconds".format(rank,endtime-starttime))
    
    if rank == 0:
        print("Pickling results into new_CMS_likes.pkl")
        with open("new_CMS_likes_finished_TESTING.pkl", 'wb') as pkl_file: 
            pickle.dump(final_loglikes,pkl_file)

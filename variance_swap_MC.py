# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 22:08:59 2022

@author: XYZW
"""


def gbm_traj(S0,r,sig,T,ns,M,q=0):
    import numpy as np
    import scipy.stats as stats
    times = np.linspace(0,T,M+1)
    arr = np.array([0]*(len(times)-1))
    cov = np.array(np.eye(len(times)-1))
    smpl = stats.multivariate_normal.rvs(arr,cov,ns)
    gbm = []
    for i in range(ns):
        traj = [S0]
        for j in range(len(times)-1):
            traj.append(traj[-1]+(r-q)*traj[-1]*(times[j+1]-times[j])+\
                        sig*traj[-1]*np.sqrt(times[j+1]-times[j])*smpl[i,j])
        gbm.append(traj)
    return np.array(gbm)

def variance_swap_price(I0,L,T,sig_K,sig,r,ns,M,q=0):
    trajectories = gbm_traj(I0,r,sig,T,ns,M,q)
    import numpy as np
    dt = T/M
    realized_vars = [np.var(np.log(trajectories[i,1::]/trajectories[i,0:-1]))*1/dt for i in range(ns)]
    payoffs = [np.exp(-r*T)*(realized_vars[i]-sig_K)*L for i in range(ns)]
    return np.mean(payoffs)


#%%

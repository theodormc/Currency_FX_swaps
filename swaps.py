# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 13:29:37 2022

@author: XYZW
"""

import sys
sys.path.append(r'C:\Users\XYZW\Documents\Python Scripts\equity exotics')

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

def currency_swap_value(FV_for,FV_dom,c_for,c_dom,r_for,r_dom,T,X0,freq = 1):
    "Principals are exchanged"
    "A domestic coupon rate is being exchanged for a foreign coupon rate"
    "Only continuous compounding is being taken into consideration"
    n = int(T*freq)
    import numpy as np
    times = np.array([1/freq * (i+1) for i in range(n)])
    pmts_for = [(c_for/freq + (i==n)) * np.exp(-times[i-1]*r_for)*FV_for 
                for i in range(1,n+1)]
    pmts_dom = [(c_dom/freq+(i==n))*np.exp(-times[i-1]*r_dom)*FV_dom 
                for i in range(1,n+1)]
    net_CFs = 1/X0*np.array(pmts_for)-pmts_dom
    price_swap = sum(net_CFs)
    return price_swap

def FX_swap_value(FV_for,FV_dom,c_for,c_dom,r_for,r_dom,T,X0,freq = 1):
    "FX swap = no principal is being exchanged"
    n = int(T*freq)
    import numpy as np
    times = np.array([1/freq * (i+1) for i in range(n)])
    pmts_for = [c_for/freq * np.exp(-times[i-1]*r_for)*FV_for 
                                        for i in range(1,n+1)]
    pmts_dom = [c_dom/freq * np.exp(-times[i-1]*r_dom)*FV_dom 
                                        for i in range(1,n+1)]
    net_CFs = 1/X0*np.array(pmts_for)-pmts_dom
    price_swap = sum(net_CFs)
    return price_swap

def variance_swap_price(I0,L,T,sig_K,sig,r,ns,M,q=0):
    trajectories = gbm_traj(I0,r,sig,T,ns,M,q)
    import numpy as np
    dt = T/M
    realized_vars = [np.var(np.log(trajectories[i,1::]/trajectories[i,0:-1]))*1/dt for i in range(ns)]
    payoffs = [np.exp(-r*T)*(realized_vars[i]-sig_K)*L for i in range(ns)]
    return np.mean(payoffs)

def equity_swap_value_MC(S1,S2,r,sig1,sig2,ro,T,ns):
    from scipy.stats import multivariate_normal as mvn
    import numpy as np
    arr = np.array([0,0])
    cov = np.array([[1,ro],[ro,1]])
    smpl = mvn.rvs(arr,cov,ns)
    aux = np.exp(np.sqrt(T)*np.array([smpl[:,0]*sig1,smpl[:,1]*sig2]).T)
    terminal_value = [S1*np.exp((r-sig1**2/2)*T),S2*np.exp((r-sig2**2/2)*T)]*aux
    payoffs = (terminal_value[:,1]-terminal_value[:,0])*(terminal_value[:,1]>terminal_value[:,0])
    price = np.exp(-r*T)*np.mean(payoffs)
    deltas1 = [aux[i,0] * (payoffs[i]>0) for i in range(ns)]
    deltas2 = [-aux[i,1] * (payoffs[i]>0) for i in range(ns)]
    delta1 = np.mean(deltas1)*np.exp(-r*T)
    delta2 = np.mean(deltas2)*np.exp(-r*T)
    return price,delta1,delta2

def maximum_equity_swap_value(S1,S2,notional,r,sig1,sig2,T,ns,M,ro):
    from bivariate_barrier_option import gbm_traj_bivariate
    gbm1,gbm2 = gbm_traj_bivariate(S1,S2,r,sig1,sig2,T,ns,M,ro)
    import numpy as np
    gbm1,gbm2 = np.array(gbm1),np.array(gbm2)
    returns1 = gbm1[:,1:M]/gbm1[:,0:M-1]-1
    returns2 = gbm2[:,1:M]/gbm2[:,0:M-1]-1
    payoffs = [(max(returns1[i,:])-max(returns2[i,:]))*notional for i in range(ns)]
    return np.mean(payoffs)*np.exp(-r*T)

def equity_swap_value(S1,S2,r,q1,q2,sig1,sig2,ro,T):
    import numpy as np
    from scipy.stats import norm
    sig_adj = np.sqrt(sig1**2+sig2**2-2*ro*sig1*sig2)
    d1 = (np.log(S2/S1)+(q1-q2-sig_adj**2/2)*T)/(sig_adj*np.sqrt(T))
    d2 = d1-np.sqrt(T)*sig_adj
    value = S2*np.exp(-q2*T)*norm.cdf(d1)-S1*np.exp(-q1*T)*norm.cdf(d2)
    print(d1,d2)
    return value

def currency_swap_floating_floating(FV_dom,FV_for,rates_dom,rates_for,\
                                    terms_dom,terms_for,X0,freq,T):
    """
    INPUTS:
        
    FV_dom = face value in domestic currency (to be paid)
    
    FV_for = face value in foreign currency (to be paid)
    """
    import numpy as np
    n = int(T*freq)
    times = [i/freq for i in range(1,n+1)]
    f1 = lambda x:np.interp(x,terms_dom,rates_dom)
    f2 = lambda x:np.interp(x,terms_for,rates_for)
    applied_rates_dom  = [f1(x) for x in times]
    applied_rates_for =  [f2(x) for x in times]
    forward_rates_dom = [applied_rates_dom[0]]+[(applied_rates_dom[i+1]*times[i+1]-\
                         applied_rates_dom[i]*times[i])/(1/freq) for i in range(len(times)-1)]
    forward_rates_for = [applied_rates_for[0]]+[(applied_rates_for[i+1]*times[i+1]-\
                         applied_rates_for[i]*times[i])/(1/freq) for i in range(len(times)-1)]
    CF_dom = [(i==n)*1/freq + forward_rates_dom[i] for i in range(len(times))]
    CF_for = [(i==n)*1/freq+forward_rates_for[i] for i in range(len(times))]
    forward_exch_rates = X0*np.exp((np.array(applied_rates_dom)-np.array(applied_rates_for))*\
                                   np.array(times))
    dom_leg = sum([CF_dom[i]*forward_exch_rates[i]*np.exp(-applied_rates_for[i]*times[i]) \
                   for i in range(len(times))])*FV_dom
    for_leg = sum([CF_for[i]*np.exp(-applied_rates_for[i]*times[i]) \
                   for i in range(len(times))])*FV_for
    value = dom_leg-for_leg
    return value
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 08:56:15 2022

@author: XYZW
"""
import sys
sys.path.append(r'C:\Users\XYZW\Documents\Python Scripts\equity exotics')
from scipy.stats import multivariate_normal as mvn
import numpy as np
import multivariate_gbm as mv_gbm
import option_price_BS as opt_BS
import pandas as pd
import sklearn.linear_model as LM
#%%
def american_equity_swap(S1,S2,r,sig1,sig2,T,ns,ro,q1 = 0,q2 = 0):
    mvar_gbm = mv_gbm.multivariate_gbm([S1,S2],S1,r,T,[sig1,sig2],ns,M,ro,[0,0])
    
#%%

#%%
def equity_swap_value_MC(S,mvn_sample,r,T,chg = 1):
    sigs = np.std(mvn_sample,axis = 0)
    aux = np.exp(np.sqrt(T)*np.array(mvn_sample*sigs))
    terminal_vals = S*np.exp((r-sigs**2/2)*T)*aux
    payoffs = (terminal_vals[:,1]-terminal_vals[:,0]*chg)*(terminal_vals[:,1]>
              terminal_vals[:,0]*chg)
    price = np.exp(-r*T)*np.mean(payoffs)
    ns = np.shape(mvn_sample)[0]
    delta1 = np.mean([aux[i,0] * (payoffs[i]>0) for i in range(ns)])*np.exp(-r*T)
    delta2 = np.mean([-aux[i,1] * (payoffs[i]>0) for i in range(ns)])*np.exp(-r*T)
    return price,delta1,delta2
#%%
def terminal_value_MC(S1,S2,r,sig1,sig2,ro,T,ns):
    arr = np.array([0,0])
    cov = np.array([[1,ro],[ro,1]])
    smpl = mvn.rvs(arr,cov,ns)
    aux = np.exp(np.sqrt(T)*np.array([smpl[:,0]*sig1,smpl[:,1]*sig2]).T)
    terminal_vals = [S1*np.exp((r-sig1**2/2)*T),S2*np.exp((r-sig2**2/2)*T)]*aux
    return terminal_vals

def equity_swap_value(S1,S2,r,sig1,sig2,ro,T,ns,chg = 1):
    
    arr = np.array([0,0])
    cov = np.array([[1,ro],[ro,1]])
    smpl = mvn.rvs(arr,cov,ns)
    aux = np.exp(np.sqrt(T)*np.array([smpl[:,0]*sig1,smpl[:,1]*sig2]).T)
    terminal_value = [S1*np.exp((r-sig1**2/2)*T),S2*np.exp((r-sig2**2/2)*T)]*aux
    payoffs = (terminal_value[:,1]-terminal_value[:,0]*chg)*(terminal_value[:,1]>
              terminal_value[:,0]*chg)
    price = np.exp(-r*T)*np.mean(payoffs)
    deltas1 = [aux[i,0] * (payoffs[i]>0) for i in range(ns)]
    deltas2 = [-aux[i,1] * (payoffs[i]>0) for i in range(ns)]
    delta1 = np.mean(deltas1)*np.exp(-r*T)
    delta2 = np.mean(deltas2)*np.exp(-r*T)
    return price,delta1,delta2

#%%
def maximum_equity_swap_value(S1,S2,notional,r,sig1,sig2,T,ns,M,ro):
    
    multivar_gbm = mv_gbm.multivariate_gbm([S1,S2],S1,r,T,[sig1,sig2],ns,M,ro,[0,0])
    returns1 = multivar_gbm[0,:,1:M]/multivar_gbm[0,:,0:M-1]-1
    
    returns2 = multivar_gbm[1,:,1:M]/gbm2[1,:,0:M-1]-1
    payoffs = [(max(returns1[i,:])-max(returns2[i,:]))*notional for i in range(ns)]
    return np.mean(payoffs)*np.exp(-r*T)
        
def variance_swap(S1,S2,r,sig1,sig2,T,ns,M,ro):
    "Exchange the variance of two equities"
    multivar_gbm = mv_gbm.multivariate_gbm([S1,S2],S1,r,T,[sig1,sig2],ns,M,ro,[0,0])
    returns1 = multivar_gbm[0,:,1:M]/multivar_gbm[0,:,0:M-1]-1
    
    returns2 = multivar_gbm[1,:,1:M]/gbm2[1,:,0:M-1]-1
    realized_variances = [(np.var(returns1[i,:])*M-np.var(returns2[i,:])*M) 
                          for i in range(ns)]
    return np.mean(realized_variances)*np.exp(-r*T)*(S1+S2)


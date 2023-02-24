# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 17:50:34 2022

@author: XYZW
"""

import numpy as np
def currency_swap_value(FV_for,FV_dom,c_for,c_dom,r_for,r_dom,T,X0,freq = 1):
    """Principals are exchanged
    A domestic coupon rate is being exchanged for a foreign coupon rate
    Only continuous compounding is being taken into consideration"""
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
    """
    FX swap = no principal is being exchanged
    
    
    """
    
    n = int(T*freq)
    import numpy as np
    times = np.array([1/freq * (i+1) for i in range(n)])
    pmts_for = [c_for/freq * np.exp(-times[i-1]*r_for)*FV_for for i in range(1,n+1)]
    pmts_dom = [c_dom/freq * np.exp(-times[i-1]*r_dom)*FV_dom for i in range(1,n+1)]
    net_CFs = 1/X0*np.array(pmts_for)-pmts_dom
    price_swap = sum(net_CFs)
    return price_swap
    

#%%
def currency_swap_value2(FV_for,FV_dom,c_for,c_dom,r_for,r_dom,T,X0,freq = 1,comp = 0):
    """Compared to the first function, this function takes into account 
    the frequency of compounding"
    
    INPUTS:
        
        X0 = no.of units of domestic currency/foreign
        
        FV_for = 
        """
        
        
    n = int(T*freq)
    times = [T-i*1/freq for i in range(n+1)][::-1] if n!=T*freq \
            else [T-i*1/freq for i in range(n)][::-1]
    import numpy as np
    if comp == 0:
        pmts_dom  =[(c_dom/freq + (i==len(times)-1))*FV_dom*np.exp(-r_dom*times[i]) for i in range(len(times))]
        pmts_for = [(c_for/freq + (i==len(times)-1))*FV_for*np.exp(-r_for*times[i]) for i in range(len(times))]
    else:
        pmts_dom  =[(c_dom/freq + (i==len(times)-1))*FV_dom*1/(1+r_dom)**times[i] for i in range(len(times))]
        pmts_for = [(c_for/freq + (i==len(times)-1))*FV_for*1/(1+r_for)**times[i] for i in range(len(times))]
    price_swap = sum(pmts_for)*1/X0-sum(pmts_dom)
    print(sum(pmts_for)*1/X0,sum(pmts_dom))
    return price_swap



#%%
def inflation_bond(FV,coupon,y,freq,CPIs,expiry):
    """
    Price of an inflation-linked bond when interest rate is linked to inflation
    
    Parameters
    ---------
    FV : face value
    
    coupon: coupon rate
    
    y: yield
    
    freq: frequency of payments
    
    CPIs: must be of length equal to the length of cash-flow payments
    
    expiry: time to maturity
    
    """
    N = int(expiry*freq)
    times = [(i+1)/freq for i in range(N)]
    interests = [coupon/freq*(1+CPIs[i]/CPIs[0])*FV for i in range(1,len(CPIs))]
    cash_flows = [interests[i]+FV*(i==len(interests)-1) for i in range(len(interests))]
    NPV = sum([np.exp(-y*times[i])*cash_flows[i] for i in range(len(cash_flows))])
    return NPV,cash_flows

def bond_test():
    FV,coupon,y,freq,expiry = 1000,0.06,0.037,4,3
    CPIs = [1000,1050,1100,1150,1200,1250,1300,1350,1400,1375,1390,1360,1390]
    print(inflation_bond(FV,coupon,y,freq,CPIs,expiry))

bond_test()
#%%

    
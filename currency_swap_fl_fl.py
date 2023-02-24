# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 09:35:08 2022

@author: XYZW
"""


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
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 18:09:49 2022

@author: XYZW
"""


def basis_swap_test():
    import numpy as np
    SOFR_rates,LIBOR_rates = [0.03,0.035,0.04],[0.022,0.027,0.034]
    SOFR_terms,LIBOR_terms = [0.5,1,1.5],[0.5,1,1.5]
    SOFR_forward_rates = [SOFR_rates[0] if i==0 else (SOFR_rates[i]*SOFR_terms[i]-\
            SOFR_rates[i-1]*SOFR_terms[i-1])/(SOFR_terms[i]-SOFR_terms[i-1])\
            for i in range(len(SOFR_rates))]
    LIBOR_forward_rates = [LIBOR_rates[0] if i==0 else (LIBOR_rates[i]*LIBOR_terms[i]-\
            LIBOR_rates[i-1]*LIBOR_terms[i-1])/(LIBOR_terms[i]-LIBOR_terms[i-1])\
            for i in range(len(SOFR_rates))]
    print(LIBOR_forward_rates)
    sa_comp_SOFR = [2*(np.exp(SOFR_forward_rates[i]/2)-1) for i in range(len(LIBOR_terms))]
    sa_comp_LIBOR = [2*(np.exp(LIBOR_forward_rates[i]/2)-1) for i in range(len(LIBOR_terms))]
    print(sa_comp_SOFR)
    print(sa_comp_LIBOR)
    notional = 1000000
    cash_flows_SOFR = [sa_comp_SOFR[i]/2*notional for i in range(len(sa_comp_SOFR))]
    cash_flows_LIBOR = [sa_comp_LIBOR[i]/2*notional for i in range(len(sa_comp_LIBOR))]
    SOFR_leg = np.dot(cash_flows_SOFR,np.exp(-np.multiply(SOFR_rates,SOFR_terms)))
    LIBOR_leg = np.dot(cash_flows_LIBOR,np.exp(-np.multiply(LIBOR_rates,LIBOR_terms)))
    print(SOFR_leg)
    print(LIBOR_leg)
    "RECEIVE SOFR PAY LIBOR"
    print(SOFR_leg-LIBOR_leg)
basis_swap_test()
#%%
def currency_swap_test():
    """
    
    """
    import numpy as np
    FV_usd,FV_eur = 5*10**6,4.5*10**6
    coupon_eur, coupon_usd = 0.03,0.02
    rate_SOFR,term_SOFR = [0.0175,0.025,0.036],[1,2,3]
    rate_EONIA,term_EONIA = [0.008,0.011,0.019],[1,2,3]
    X0 = 1.2
    forward_rates = X0*np.exp((np.array(rate_EONIA)-np.array(rate_SOFR)*np.array(term_SOFR)))
    CF_euros = [(coupon_eur+(i==3))*FV_eur for i in range(1,4)]
    CF_dollars = [(coupon_usd+(i==3))*FV_usd for i in range(1,4)]
    equivalent_CF_euros = [CF_euros[i]*forward_rates[i] for i in range(len(CF_euros))]
    eur_leg = sum([equivalent_CF_euros[i]*np.exp(-rate_SOFR[i]*term_SOFR[i]) for i in \
                   range(len(term_SOFR))])
    usd_leg = sum([CF_dollars[i]*np.exp(-rate_SOFR[i]*term_SOFR[i]) for i in \
                   range(len(term_SOFR))])
    currency_swap_value = eur_leg-usd_leg
    print(currency_swap_value)
currency_swap_test()
#%%
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

def test_currency_swap_fl_fl():
    FV_usd,FV_eur = 5*10**6,4.5*10**6
    SOFR_rates,SOFR_terms = [0.0175,0.025,0.036],[1,2,3]
    EONIA_rates,EONIA_terms = [0.008,0.011,0.019],[1,2,3]
    X0,freq,T = 1.2,1,3
    print(currency_swap_floating_floating(FV_eur,FV_usd,EONIA_rates,SOFR_rates,\
                            EONIA_terms,SOFR_terms,X0,freq,T))
test_currency_swap_fl_fl() 
#%%
def currency_swap_floating_floating_test():
    """
    We value a floating for floating currency swap where principals are also exchanged
    The example below regards a 3-year currency swap where floating rates in USD (SOFR)
    are exchanged for floating rates in euros (EONIA)
    """
    import numpy as np
    FV_usd,FV_eur = 5*10**6,4.5*10**6
    coupon_eur, coupon_usd = 0.03,0.02
    SOFR_rates,SOFR_terms = [0.0175,0.025,0.036],[1,2,3]
    EONIA_rates,EONIA_terms = [0.008,0.011,0.019],[1,2,3]
    X0 = 1.2
    forward_exch_rates = X0*np.exp((np.array(EONIA_rates)- \
                                    np.array(SOFR_rates)*np.array(SOFR_terms)))
    forward_usd_rates = [SOFR_rates[0] if i==0 else (SOFR_rates[i]*SOFR_terms[i]-\
            SOFR_rates[i-1]*SOFR_terms[i-1])/(SOFR_terms[i]-SOFR_terms[i-1])\
            for i in range(len(SOFR_rates))]
    forward_eur_rates = [EONIA_rates[0] if i==0 else (EONIA_rates[i]*EONIA_terms[i] - \
                         EONIA_rates[i-1]*EONIA_terms[i-1])/(EONIA_terms[i]-EONIA_terms[i-1])\
            for i in range(len(EONIA_rates))]
    CF_euros = [(coupon_eur+forward_eur_rates[i-1]+(i==3))*FV_eur for i in range(1,4)]
    CF_dollars = [(coupon_usd+forward_usd_rates[i-1]+(i==3))*FV_usd for i in range(1,4)]
    equivalent_CF_euros = [CF_euros[i]*forward_exch_rates[i] for i in range(len(CF_euros))]
    eur_leg = sum([equivalent_CF_euros[i]*np.exp(-SOFR_rates[i]*SOFR_terms[i]) for i in \
                   range(len(SOFR_rates))])
    usd_leg = sum([CF_dollars[i]*np.exp(-SOFR_rates[i]*SOFR_terms[i]) for i in \
                   range(len(SOFR_terms))])
    currency_swap_value = eur_leg-usd_leg
    print(currency_swap_value, eur_leg, usd_leg)
    print(CF_euros)
    print(CF_dollars)
currency_swap_floating_floating_test()
#%%
def timing_adjustment(ro,T1,T2,vol1,vol2,rf,freq):
    import numpy as np
    adj = np.exp(-ro*(T2-T1)*vol1*vol2*rf/(1+rf/freq))
    fwd_applied_rate = rf+adj
    return fwd_applied_rate

def test_timing():
    import numpy as np
    ro,T1,T2,vol1,vol2,rf,freq = 1,9,8,0.2,0.2,0.06,1
    adj = timing_adjustment(ro,T1,T2,vol1,vol2,rf,freq)
    FV = 1000
    value = FV*(0.06+adj)*0.5*np.exp(-rf*T1)
    print(value)
    
test_timing()
#%%
def test_interp():
    import numpy as np
    times = [1,2,3,4]
    rates = [0.025,0.0375,0.045,0.0535]
    f = lambda x:np.interp(x,times,rates)
    print(f(4.5))
test_interp()
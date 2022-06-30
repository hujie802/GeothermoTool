# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:18:23 2019

@author: hujie
"""
import numpy as np
#paleo sealevel correction
def paleosealevel(lays,rho_m,rho_w):
    paleo_s = []
    time = []
    for i in range(len(lays)):
        paleo_s.append(lays[i].paleosealevel)
        time.append(lays[i].age_end)
    paleos = np.array(paleo_s)*rho_w/(rho_m-rho_w)
    return paleos,time
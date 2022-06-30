# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:14:27 2019

@author: hujie
"""
import numpy as np
def sedi_load(nlays,rho_m,rho_w):

    """
    sediment load correction
    rho_m: density of mantle
    rho_w: density of water
    """
    zz = []
    time = []
    for j in range(len(nlays)):
        tt = 0
        lays = nlays[j]
        for i in range(len(lays)):
            tt = tt+lays[i].density()*lays[i].thickness()
        rho_s = tt/lays[-1].bottom
        kk = lays[-1].bottom*(rho_m-rho_s)/(rho_m-rho_w)
        zz.append(kk)
        time.append(lays[0].age_end)
        qq = np.array(zz)
    time.append(lays[0].age_start)
    return(qq,time)

# -*- coding: utf-8 -*-
"""
Created on Mon May 27 21:21:11 2019

@author: hujie
"""
import numpy as np 
#paleopathmetry correction
def paleobathymetry(layers):
    paleobb = []
    for i in range(len(layers)):
        paleobb.append(layers[i].paleobathymetry)
    paleob = np.array(paleobb)
    paleocc = []
    for i in range(len(layers)):
        paleocc.append(layers[i].paleosealevel)
    paleoc = np.array(paleocc)
    return paleob-paleoc
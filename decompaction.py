# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 18:58:20 2019

@author: hujie
"""

import numpy as np
from numpy import exp
import copy
from Layer import Layer


# This code is for decompact the strata.
# y1-top, y2-bottom before decompation
# y3-top, y4-bottom after decompation
# phi0 - porosity in ground
# c - 1/y0  The typical values of the depth constant y 0 for common lithologies
# in the North Sea were given by Sclater and Christie (1980), who expressed
# Allen, 2013

def decomp(lays):
    nlayers = [lays,]
    layers_use = copy.deepcopy(lays)
    for j in range(len(lays)-1):
        layers_use.pop(0)
        temp_start = 0
        layers_temp = []
        for i in range(len(layers_use)): 
            layer_temp = copy.deepcopy(layers_use[i])
            y2 = layer_temp.bottom
            y1 = layer_temp.top
            y3 = temp_start 
            phi0 = layer_temp.phi0()
            c = layer_temp.c()
            #binary search
            left = 5
            right = 100000
            eps = 1 #in m
            mid = (left+right)/2
            while abs(y2-y1-phi0/c*(exp(-c*y1)-exp(-c*y2))+y3-mid+phi0/c*(exp(-c*y3)-exp(-c*mid)))>eps:
                mid = (left+right)/2
                if (y2-y1-phi0/c*(exp(-c*y1)-exp(-c*y2))+y3-left+phi0/c*(exp(-c*y3)-exp(-c*left)))*(y2-y1-phi0/c*(exp(-c*y1)-exp(-c*y2))+y3-mid+phi0/c*(exp(-c*y3)-exp(-c*mid)))<0:
                    right = mid
                    # print(left)
                else:
                    left = mid
            # print(mid)
            temp_bottom = mid
            layer_temp.top = temp_start
            layer_temp.bottom = temp_bottom
            layers_temp.append(layer_temp)
            temp_start = copy.deepcopy(temp_bottom)
        layers_use = copy.deepcopy(layers_temp)
        nlayers.append(layers_temp)
    return nlayers
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 18:51:54 2019

@author: hujie
"""
import pandas as pd 
from Layer import Layer

def read_file(filename,phi0_sand,phi0_mud,phi0_carbon,phi0_coal,c_sand,c_mud,c_carbon,c_coal,rho_sand,rho_mud,rho_carbon,rho_coal,rho_water):
    """
    the function to read files
    top - top depth of all layers
    age_end - deposit end age of all layers
    sand,mud,carbon,coal: proportion of those compositions, 0-1
    phi0_sand,phi0_mud,phi0_carbon,phi0_coal: initial porosity , 0-1
    c_sand,c_mud,c_carbon,c_coal: Compaction coefficient, m^-1
    rho_sand,rho_mud,rho_carbon,rho_coal,rho_water: density,kg/m3

    """    
    f = pd.read_csv(filename)
    top = list(f.Depth)
    top.insert(0,0)   
    top.pop()
    f['Top'] = top
    age_end = list(f.Age)
    age_end.insert(0,0)
    age_end.pop()
    f['Age_end'] = age_end
    layers = list(f.Formation)
    for i in range(len(f.Age)):
        layers[i] = Layer(f.Formation[i],f.Age[i],f.Age_end[i],f.Top[i],f.Depth[i],f.Paleobathymetry[i],f.Paleosealevel[i],
              f.Sand[i],f.Mud[i],f.Carbon[i],f.Coal[i],
              phi0_sand,phi0_mud,phi0_carbon,phi0_coal,
              c_sand,c_mud,c_carbon,c_coal,
              rho_sand,rho_mud,rho_carbon,rho_coal,rho_water)
    return layers
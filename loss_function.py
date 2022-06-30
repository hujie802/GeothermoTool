# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 19:36:41 2019

@author: hujie
"""

from T_fixed_stp_fixed_Jit_node_solver import temp_solver
import numpy as np
from scipy import interpolate
import time

def loss(G,op_points_age,age,subsidence,depth,rock,T0,Tm,density_s,density,conductivity,heat_production,heat_capacity,alpha):
    """
    loss_function: the residuals between observed tectonic subsidence and predicted tectonic subsidence 

    G: vector, strain rate,s^-1, the optimizational value; 
    op_points_age: vector,age of strain rate value, Ma;
    age: vector, age of the formations,Ma;
    subsidence:vector,observed subsidence, m;
    depth: vector,structrure of lithosphere, the depths of the upper crust, lower crust, lithosphere,asthenosphere,m
    rock: rock type of the lithosphere,0,1,2,3
    T0: scalar, surface temperature, C
    Tm: scalar, lithosphere bottom temperature, usually 1330C
    density_s: scalar,density of sediments, kg/m3
    density,conductivity,heat_production,heat_capacity,alpha : Thermophysical parameters of the lithosphere

    """
    age_end = age[len(age)-1] # the time of the initial rifting
    #interpolate G to a function
    G_func = interpolate.interp1d(op_points_age,G,kind ='linear')

    #solve the temperature function
    T,age_new,beta,HF_need,uc_depth,lc_depth,lm_depth,uc_density,lc_density,lm_density,as_density= temp_solver(T0,Tm,rock,density,conductivity,heat_production,heat_capacity,depth,alpha,G_func,age_end)
    #solve the tectonic subsidence  
    Y = (uc_density*uc_depth+lc_density*(lc_depth-uc_depth)+lm_density*(lm_depth-lc_depth)+as_density*(lm_depth[0]-lm_depth)
    -uc_density[0]*uc_depth[0]-lc_density[0]*(lc_depth[0]-uc_depth[0])-lm_density[0]*(lm_depth[0]-lc_depth[0]))/(as_density-density_s)
#    #拟合所有点
    sub_func = interpolate.interp1d(age,subsidence,kind ='linear')
#    calc the subsidence in ages
    sub_cal = sub_func(age_new)
    #use L2norm
    misfit = np.linalg.norm((Y-sub_cal),ord=2)
    return misfit
    
    
    
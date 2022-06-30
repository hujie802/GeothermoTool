# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 14:17:21 2021

@author: HP
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:30:23 2019

@author: hujie
"""

#from T_moving_node_solver import temp_solver
from T_fixed_stp_fixed_Jit_node_solver import temp_solver
import pandas as pd
import numpy as np
from scipy import optimize,interpolate
import copy
import matplotlib.pyplot as plt
from loss_function import loss
from scipy.optimize import minimize
import time
import os


def subsub(G,op_points_age,sub,depth,rock,T0,Tm,density_s,density,conductivity,heat_production,heat_capacity,alpha):
        # calc the subsidence with strain rate
        age = sub.age
        age_end = age[len(age)-1]
        #interpolate G to a function
        G_func = interpolate.interp1d(op_points_age,G,kind ='linear')
        #solve temperature function
        T,age_new,beta,HF_need,uc_depth,lc_depth,lm_depth,uc_density,lc_density,lm_density,as_density= temp_solver(T0,Tm,rock,density,conductivity,heat_production,heat_capacity,depth,alpha,G_func,age_end)

        Y = (uc_density*uc_depth[0]/beta+lc_density*(lc_depth[0]-uc_depth[0])/beta+lm_density*(lm_depth[0]-lc_depth[0])/beta+as_density*(lm_depth[0]-lm_depth[0]/beta)-uc_density[0]*uc_depth[0]
    -lc_density[0]*(lc_depth[0]-uc_depth[0])-lm_density[0]*(lm_depth[0]-lc_depth[0]))/(as_density-density_s)
        return Y,T,age_new,beta,HF_need,uc_depth,lc_depth,lm_depth,uc_density,lc_density,lm_density

def smooth(a,WSZ):
    # a:raw data，NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))


def thermo_ev(sub,T0,Tm,depth,density,conductivity,heat_production,heat_capacity,alpha,density_s):
    rock = [0,1,2,3]    #upper crust，lower crust，lithosphere,Asthenosphere
    # sub = pd.read_csv(name)
    age = sub.age
    subsidence = sub.subsidence #observed subsidenc

    # n_op = len(age)
    # op_points_age = age

    # number of optimization point
    n_op = 11
    # age array of the strain rate value
    op_points_age = np.linspace(age[len(age)-1],0,n_op) 
        
    # G0 = np.linspace(-17,-17,n_op)
    G0 = np.random.uniform(-15,-18,n_op)
    lb = -23
    ub = -14
    #bounds of the strain rate
    lw = [lb]*n_op
    up = [ub]*n_op
    bounds = list(zip(lw,up))
        
    ##    #optimization
    start_time = time.time()
    res = minimize(loss,G0,(op_points_age,age,subsidence,depth,rock,T0,Tm,density_s,density,conductivity,heat_production,heat_capacity,alpha),method = 'SLSQP',bounds = bounds,options={'maxiter': 200})
    end_time = time.time()
    print(end_time-start_time)

    #the optimized results
    strain_rate = res.x
    #calc the subsidence using the optimized strain rate
    Y,T,age_new,beta,HF_need,uc_depth,lc_depth,lm_depth,uc_density,lc_density,lm_density = subsub(strain_rate,op_points_age,sub,depth,rock,T0,Tm,density_s,density,conductivity,heat_production,heat_capacity,alpha)
    HF_for_use =  smooth(HF_need*1000,21) # smooth the heat flow
    return Y,T,age_new,beta,HF_for_use,strain_rate,op_points_age





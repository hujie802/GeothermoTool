# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 16:31:18 2019

@author: hujie
"""

import numpy as np
from numba import jit

#number of the time step
number_of_time_step = 200
#compile the function
@jit(nopython=True)
def state_thermal(D_uc0,D_lc0,D_lm0,nzz,nd,T0,Tm,conductivity,heat_production,stp,zm,NM):
    #the calculation grid is limited at lithosphere                 
    #use steady-state function to calc the initial temperature
    #
    K_temp = np.zeros(nd)
    A_temp = np.zeros(nd)
    T_m_0 = np.zeros(NM) 
    #giving value,thermal conductivity，heat production
    K_temp[nzz <D_uc0] = conductivity[0]
    K_temp[((nzz >= D_uc0) & (nzz <D_lc0))] = conductivity[1]
    K_temp[nzz >= D_lc0] = conductivity[2]
    
    A_temp[nzz < D_uc0] = heat_production[0]
    A_temp[((nzz >= D_uc0) & (nzz < D_lc0))] = heat_production[1]
    A_temp[nzz >= D_lc0] = heat_production[2]
    
    #Function left,right
    L_s = np.zeros((nd,nd))
    R_s = np.zeros(nd)
    for i in range(nd):
        #upper boundary
        if i==0:
            L_s[i,i] = 1
            R_s[i]=T0
        #bottom boundary
        elif i==nd-1:
            L_s[i,i] = 1
            R_s[i]=Tm
        else:
            L_s[i,i-1] = -(K_temp[i]+K_temp[i-1])/2/stp**2
            L_s[i,i+1] = -(K_temp[i+1]+K_temp[i])/2/stp**2
            L_s[i,i] = (K_temp[i+1]+K_temp[i])/2/stp**2+(K_temp[i]+K_temp[i-1])/2/stp**2            
            R_s[i] = A_temp[i]
            
    #solve functions
    T_start = np.linalg.solve(L_s,R_s)
    HF0 = (T_start[3]-T_start[1])/stp/2*conductivity[0]+A_temp[0]*stp
    
    #node to marker
    for z in range(NM):
        y = int(zm[z]/stp)      
        q = zm[z]/stp-y    # weighting value
        T_m_0[z] = T_m_0[z]+ T_start[y]*(1-q)+T_start[y+1]*q
        
    return T_start,HF0,T_m_0

#compile the function
@jit(nopython=True)
def transient_therm(dt,nd,NM,zm,Rock_m,Rho_m,K_m,T_m,A_m,Cp_m,Alpha_m,stp,T0,Tm):
#initialize the grid
    Rock = np.zeros(nd)      #rock
    Rho = np.zeros(nd)       #density
    K = np.zeros(nd)         #thermal conductivity 
    T = np.zeros(nd)         #Temperature
    A = np.zeros(nd)         #Heat production
    Cp = np.zeros(nd)        #Specific capacity
    Alpha = np.zeros(nd)     #Coefficient of thermal expansion
    Weight = np.zeros(nd)    #weighting value
           
    #marker to node      
    for im in range(NM):
        i = int(zm[im]/stp)  #marker 
        # calc the distance between marker to node
        dz = zm[im]/stp-i
        #-+--------[xn]-----(im)----+--------[xn+1]--------+--------
        if(dz<0.5):
            Rock[i] = Rock[i]+Rock_m[im]*(1-dz) 
            Rho[i] = Rho[i]+Rho_m[im]*(1-dz)
            K[i] = K[i]+K_m[im]*(1-dz)
            T[i] = T[i]+T_m[im]*(1-dz)
            A[i] = A[i]+A_m[im]*(1-dz)
            Cp[i] = Cp[i]+Cp_m[im]*(1-dz)
            Alpha[i] = Alpha[i]+Alpha_m[im]*(1-dz)
            Weight[i] = Weight[i]+(1-dz)
            
        else:
            Rock[i+1] = Rock[i+1]+Rock_m[im]*dz 
            Rho[i+1] = Rho[i+1]+Rho_m[im]*dz
            K[i+1] = K[i+1]+K_m[im]*dz
            T[i+1] = T[i+1]+T_m[im]*(dz)
            A[i+1] = A[i+1]+A_m[im]*(dz)
            Cp[i+1] = Cp[i+1]+Cp_m[im]*dz
            Alpha[i+1] = Alpha[i+1]+Alpha_m[im]*dz
            Weight[i+1] = Weight[i+1]+dz
#   calc weighting value       
    Rho = Rho/Weight
    K = K/Weight
    T = T/Weight
    A = A/Weight
    Alpha = Alpha/Weight
    Cp = Cp/Weight
    Rock = Rock/Weight
    #set upper point
    T[0] = T0
#    T[-1] = Tm
#left part: Save the temperature after moving
    L = np.zeros((nd,nd))
# right part 
    R = np.zeros(nd)
    for h in range(nd):
#      upper boundary
        if h == 0:
            L[0,0] = 1
            R[0] = T0
        #bottom boundary
        elif h == nd-1:
            L[nd-1,nd-1] = 1
            R[nd-1]= T[-1]
#       #internal node
        else:
            L[h,h-1] = -(K[h-1]+K[h])/stp**2/2         # T+dt[i-1]
            L[h,h] = Rho[h]*Cp[h]/dt + (K[h+1]+K[h])/stp**2/2 + (K[h-1]+K[h])/stp**2/2 #T+dt[i]
            L[h,h+1] = -(K[h+1]+K[h])/stp**2/2          # T+dt[i+1] 
#                
            R[h] = Rho[h]*Cp[h]/dt*T[h]+A[h]    
#    #  solve equation      
    T_new = np.linalg.solve(L,R)      
    T_delta = T_new-T
    #node to marker
    T_m_new = np.zeros(NM)
    for z in range(NM):
        y = int(zm[z]/stp)      #marker to former node
        q = zm[z]/stp-y    #distance to former node
        T_m_new[z] = T_m[z]+ T_delta[y]*(1-q)+T_delta[y+1]*q    
        
    return T_new,T_m_new, Rho,Rock

#@jit(nopython=False)
def temp_solver(T0,Tm,rock,density,conductivity,heat_production,heat_capacity,depth,alpha,G_func,age_end):
    #
    T0 = T0      #upper boundary
    Tm = Tm    #Lithospheric bottom boundary temperature
    
    D_uc0 = depth[0] #depth of upper crust m
    D_lc0 = depth[1] #depth of lower crust m
    D_lm0 = depth[2] #depth of lithosphere crust m

    # Solving equations using one-dimensional finite difference method
    # define the initial drid
    nd = 101  #number of node 
    stp = D_lm0/(nd-1) #step 
    nzz = np.linspace(0,D_lm0,nd) #grid
    
    ##define initial marker
    nm = 4         # 1 step include 4Marker
    dm = stp/nm     #marker step 
    NM = nm*(nd-1)  #number of marker
    zm = np.linspace(dm*0.5,D_lm0-dm*0.5,NM) #initial marker grid
    
    T_start,HF0,T_m0 = state_thermal(D_uc0,D_lc0,D_lm0,nzz,nd,T0,Tm,np.array(conductivity),np.array(heat_production),stp,zm,NM)

    #----giving value to marker ---   
    #initize the marker
    Rock_m = np.zeros(NM) #rock type
    K_m = np.zeros(NM)    #thermal conductivity
    A_m = np.zeros(NM)    #heat production
    Rho_m = np.zeros(NM)  #density
    Cp_m = np.zeros(NM)   #heat capacity
    T_m = np.zeros(NM)    #temperature
    Alpha_m = np.zeros(NM) #Coefficient of thermal expansion
    Rho_m_raw = np.zeros(NM) #raw density
    
    Rock_m[zm < D_uc0] = rock[0]
    Rock_m[((zm >=D_uc0) & (zm< D_lc0))] = rock[1]
    Rock_m[zm >=D_lc0] = rock[2]
    
    K_m[zm < D_uc0] = conductivity[0]
    K_m[((zm >= D_uc0) & (zm < D_lc0))] = conductivity[1]
    K_m[zm >=D_lc0] = conductivity[2]

    A_m[zm < D_uc0] = heat_production[0]
    A_m[((zm >=D_uc0) & (zm< D_lc0))] = heat_production[1]
    A_m[zm >=D_lc0] =heat_production[2]
    #原始密度，没考虑热膨胀    
    Rho_m_raw[zm < D_uc0] = density[0]
    Rho_m_raw[((zm >=D_uc0) & (zm< D_lc0))] = density[1]
    Rho_m_raw[zm >=D_lc0] = density[2]
    
    Cp_m[zm < D_uc0] = heat_capacity[0]
    Cp_m[((zm >=D_uc0) & (zm < D_lc0))] = heat_capacity[1]
    Cp_m[zm >=D_lc0] = heat_capacity[2]
    
    Alpha_m[zm < D_uc0] = alpha[0]
    Alpha_m[((zm >=D_uc0) & (zm < D_lc0))] = alpha[1]
    Alpha_m[zm >=D_lc0] = alpha[2]
    
    T_m = T_m0    
    #solve density
    Rho_m = Rho_m_raw*(1-Alpha_m*T_m)
        
    distance_m0 = zm[1:]-zm[:-1]  #marker step 
    Rho_m_mid0 = (Rho_m[:-1]+Rho_m[1:])/2 #marker density
    gm0 = distance_m0*Rho_m_mid0    #quality
    gm0= np.insert(gm0,0,0)
    
    uc_density0 = gm0[Rock_m==rock[0]].sum()/(zm[Rock_m==rock[0]][-1]-zm[Rock_m==rock[0]][0])
    lc_density0 = gm0[Rock_m==rock[1]].sum()/(zm[Rock_m==rock[1]][-1]-zm[Rock_m==rock[0]][-1])
    lm_density0 = gm0[Rock_m==rock[2]].sum()/(zm[Rock_m==rock[2]][-1]-zm[Rock_m==rock[1]][-1])

    as_density0 = density[3]*(1-Tm*alpha[3])
    
#    print(as_density)
    # save density      
    uc_density = [uc_density0] #upper crust
    lc_density = [lc_density0] #lower crust
    lm_density = [lm_density0] #lithosphere mantle
    as_density = [as_density0] #asothenosphere
    
#    print(uc_density0,lc_density0,lm_density0)
    
    need_T =   [T_start] #save temperature
    age_new =  [age_end] #save time
    uc_depth = [D_uc0] #save depth of upper crust
    lc_depth = [D_lc0] #save depth of lower crust
    lm_depth = [D_lm0] #save depth of lithoshpere mantle
    
    HF_need =  [HF0]    #save heat flow
    #start loop
    dt =  age_end/number_of_time_step*1e6*365*24*3600
    t = dt #initial time step

    beta = [1] # Stretching factor

    for i in range(number_of_time_step):      
#        print(i)
        
        time_temp = age_end-t/(1e6*365*24*3600) #age inversion, Ma
#        print(t/(1e6*365*24*3600))
        if int(time_temp*1e6*365*24*3600)<1:   #Allow single digit error
            time_temp = 0

        G = 10**G_func(time_temp) #strain rage
##      # move depth
        D_lm0 = D_lm0-G*D_lm0*dt  #lithosphere bottom depth
        D_lc0 = D_lc0-G*D_lc0*dt  #lower crust bottom depth
        D_uc0 = D_uc0-G*D_uc0*dt  #upper crust bottom depth
   
        #move marker
        vzm = G*zm
        zm = zm-vzm*dt
        
#        beta_org = beta_org+G*dt
        beta.append(depth[2]/D_lm0)
        #save to list
        uc_depth.append(D_uc0)
        lc_depth.append(D_lc0)
        lm_depth.append(D_lm0)
#        print(zm[-1])
        n_half_stp = int((depth[2]-zm[-1])/stp*2)
        if n_half_stp>0:
            start_add = depth[2]-n_half_stp/2*stp+stp/nm/2
            end_add = depth[2]-stp/nm/2
            n_add_marker = int(n_half_stp*nm/2)
#            print(n_add_marker)
            zm_add = np.linspace(start_add,end_add,n_add_marker)
            zm = np.append(zm,zm_add)
            Rock_m = np.append(Rock_m,[rock[3]]*(n_add_marker))
            Rho_m = np.append(Rho_m,[as_density0]*(n_add_marker))
            K_m = np.append(K_m,[conductivity[3]]*(n_add_marker))
            T_m = np.append(T_m,[Tm]*(n_add_marker))
            A_m = np.append(A_m,[heat_production[3]]*(n_add_marker))
            Cp_m = np.append(Cp_m,[heat_capacity[3]]*(n_add_marker))
            Alpha_m = np.append(Alpha_m,[alpha[3]]*(n_add_marker))
            Rho_m_raw = np.append(Rho_m_raw,[density[3]]*(n_add_marker))

        NM = len(zm) 
#        print(Rho_m)
        T_new,T_m_new ,Rho,Rock= transient_therm(dt,nd,NM,zm,Rock_m,Rho_m,K_m,T_m,A_m,Cp_m,Alpha_m,stp,T0,Tm)
        need_T.append(T_new)
        T = T_new
        T_m = T_m_new
#        print(T_new)
#        print(stp)
        HF = (T[3]-T[1])/stp/2*conductivity[0] +A_m[0]*stp #calc heat flow
        HF_need.append(HF)
        
#       density change
        Rho_m = Rho_m_raw*(1-Alpha_m*T_m)
        
        distance_m = zm[1:]-zm[:-1]  #marker step
        Rho_m_mid = (Rho_m[:-1]+Rho_m[1:])/2 #marker density
        gm = distance_m*Rho_m_mid    #quality
        gm = np.insert(gm,0,0)
#        print(distance_m)
#        Rho_func = interpolate.interp1d(zm,Rho_m,kind='linear')
        #average density of upper crust
        uc_density.append(gm[Rock_m==rock[0]].sum()/(zm[Rock_m==rock[0]][-1]-zm[Rock_m==rock[0]][0]))
#        print(uc_density)
#        #average density of lower crust
        lc_density.append(gm[Rock_m==rock[1]].sum()/(zm[Rock_m==rock[1]][-1]-zm[Rock_m==rock[0]][-1]))
#        print(lc_density)
##        #average density of lithosphere mantle
        lm_density.append(gm[Rock_m==rock[2]].sum()/(zm[Rock_m==rock[2]][-1]-zm[Rock_m==rock[1]][-1]))
#        print(uc_density)

        if (list(Rho[Rock==rock[3]])== []):
            as_density.append(as_density0)
        else:
            as_density.append(gm[Rock_m==rock[3]].sum()/(zm[Rock_m==rock[3]][-1]-zm[Rock_m==rock[2]][-1]))

        age_new.append(time_temp)
        t = t+dt
    return need_T,age_new,np.array(beta),np.array(HF_need),np.array(uc_depth),np.array(lc_depth),np.array(lm_depth),np.array(uc_density),np.array(lc_density),np.array(lm_density),np.array(as_density)

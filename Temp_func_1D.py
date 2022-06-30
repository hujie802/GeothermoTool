import numpy as np

def func_1D(depth_layer,k_layer,A_layer,T0,bottom_bound,boundary_cond):
    #  depth_layer,k_layer,A_layer :depth of each layer in m,thermal conductivity of each layer in W/mK,heat production in W/m3
    # T0 is upper boundary in C
    #bottom_bound: bottom boundary, in C or W/m2
    #boundary_cond: HF or Temp
    cal_depth_top = 0
    cal_depth_bottom = depth_layer[len(depth_layer)-1]
    #++++++++++++++++++parameters++++++++++++++++
    nd = 1001 #number of node
    stp = (cal_depth_bottom-cal_depth_top)/(nd-1) # step
    nz = np.linspace(cal_depth_top,cal_depth_bottom,nd)
    #++++++++++++++++++parameters++++++++++++++++

    #giving values
    K_temp = np.zeros(nd)
    A_temp = np.zeros(nd)
    for i,d in enumerate(depth_layer):
        if i==0:
            top_lim =0      # top limit of the slice
        else:
            top_lim = depth_layer[i-1]  # bottom limit of the slice
        bot_lim = depth_layer[i]
        K_temp[(nz <= bot_lim) & (nz >= top_lim)] = k_layer[i]  # giving value
        A_temp[(nz <= bot_lim) & (nz >= top_lim)] = A_layer[i]  # giving value

    L_s = np.zeros((nd,nd))
    R_s = np.zeros((nd))
    #

    for i in range(nd):
        if i == 0:   #upper boundary
            L_s[i,i] = 1
            R_s[i] = T0
        elif i == nd-1:  #bottom boundary
            if boundary_cond == 'HF':
            #heat flow boundary
               L_s[i,i-1] = -1*K_temp[i]/stp
               L_s[i,i] = 1*K_temp[i]/stp
               R_s[i] = bottom_bound # mW/m2 to W/m2
            else:#temperature boundary
                L_s[i,i] = 1
                R_s[i] = bottom_bound
        else:
            L_s[i,i-1] = -(K_temp[i]+K_temp[i-1])/2/stp**2
            L_s[i,i+1] = -(K_temp[i+1]+K_temp[i])/2/stp**2
            L_s[i,i] = (K_temp[i+1]+K_temp[i])/2/stp**2+(K_temp[i]+K_temp[i-1])/2/stp**2
            R_s[i] = A_temp[i]

    T_start = np.linalg.solve(L_s,R_s)
    return nz,T_start

if __name__ == "__main__":
    # +++++++++++++++++input values++++++++++++++
    depth_layer = [2000, 4000, 8000, 20000]  # depth of each layer in m
    k_layer = [2, 3, 3, 3]  # thermal conductivity of each layer in W/mK
    A_layer = [1.5e-6, 2e-6, 3e-6, 4e-6]  # heat production in W/m3

    # boundary condition
    T0 = 0  # top boundary in C
    # Tb = 500  # bottom boundary in C
    # bottom_bound = Tb
    # boundary_cond = 'TEMP'
    bottom_bound = 30
    boundary_cond = 'HF'
    nz,tep = func_1D(depth_layer, k_layer, A_layer, T0, bottom_bound, boundary_cond)
import numpy as np
import pandas as pd
from scipy import interpolate
import scipy.sparse as sparse
import matplotlib.pyplot as plt

#

def smooth(a,WSZ):
    # a:raw data，NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))

def func_2D(xlim,ylim,layer_funcs,K_layer,A_layer,upper_bounds,bottom_bounds,boundary_cond):
    #xlim,ylim the calc range of x and y in m, note the top cal boundary is 0m, the left boundry is 0m
    #  layer_funcs,k_layer,A_layer :depth func of each layer in m,thermal conductivity of each layer in W/mK,heat production in W/m3
    # T0 is upper boundary in C
    # upper_bounds: upper boundary, in C
    #bottom_bounds: bottom boundary, in C or W/m2
    #boundary_cond: HF or Temp

    #++++++++++++++++++parameters++++++++++++++++
    nx = 201
    ny = 201#number of node
    xstp = xlim/(nx-1) # xstep
    ystp = ylim/(ny-1) # ystep
    x_arr = np.linspace(0,xlim,nx)
    y_arr = np.linspace(0,ylim,ny)

    bottom_bounds = bottom_bounds(x_arr)
    #++++++++++++++++++giving values++++++++++++++++
    K_temp = np.zeros((ny,nx))
    A_temp = np.zeros((ny,nx))
    for col in range(nx):  # giving value by column
        x_tp = x_arr[col]  #giving vaule to a column
        for l,fuc in enumerate(layer_funcs):
            if l == 0:
                top_lim = 0  # top limit of the slice
            else:
                top_lim = float(layer_funcs[l-1](x_tp))  # bottom limit of the slice
            bot_lim = float(fuc(x_tp))
            K_temp[(y_arr <= bot_lim) & (y_arr >= top_lim),col] = K_layer[l]  # giving value
            A_temp[(y_arr <= bot_lim) & (y_arr >= top_lim),col] = A_layer[l]  # giving value

    L_s = sparse.dok_matrix((nx*ny,nx*ny))
    R_s=np.zeros(nx*ny)
    #
    for i in range(ny):
        for j in range(nx):
            # Composing matrix of coefficients L()
            # and vector (column) of right parts R()
            # Globe index
            k = j*ny+i
            if i == 0:   #upper boundary
                L_s[k,k] = 1
                R_s[k] = upper_bounds
            elif i == ny-1:  #bottom boundary
                if boundary_cond == 'HF':#Constant flux
                #heat flow boundary
                    L_s[k,k-1] = -K_temp[i,j]/ystp
                    L_s[k,k] = K_temp[i,j]/ystp
                    R_s[k] = bottom_bounds[j]
                else: #Constant temperature boundary
                    L_s[k,k] = 1
                    R_s[k] = bottom_bounds[j]
            elif (j == 0 and i > 0 and i < ny - 1):
            # Insulating boundary: T(i,j)=T(i,j+1)
                L_s[k, k] = 1
                L_s[k, k + ny] = -1
                R_s[k] = 0
            # right boundary
            elif (j == nx-1 and i > 0 and i < ny - 1):
            # Insulating boundary: T(i,j-1)=T(i,j)
                L_s[k, k-ny] = 1
                L_s[k, k] = -1
                R_s[k] = 0
            else:
                L_s[k, k - ny] = -(K_temp[i, j - 1] + K_temp[i, j])/2/ xstp**2
                L_s[k, k + ny] = -(K_temp[i, j + 1] + K_temp[i, j])/2/ xstp**2
                L_s[k, k - 1] = -(K_temp[i - 1, j] + K_temp[i, j])/2/ystp**2
                L_s[k, k + 1] = -(K_temp[i + 1, j] + K_temp[i, j])/2/ystp**2
                L_s[k, k] =  (K_temp[i, j] + K_temp[i, j + 1])/2/xstp**2 \
                            + (K_temp[i, j] + K_temp[i, j - 1])/2/xstp**2 \
                            + (K_temp[i, j] + K_temp[i + 1, j])/2/ystp**2 \
                            + (K_temp[i, j] + K_temp[i - 1, j])/2/ystp**2
                R_s[k] = A_temp[i, j]
    # solve the eqoution
    L_s = L_s.tocsc()  # csc format is more efficient
    from scipy.sparse.linalg import spsolve
    S = spsolve(L_s, R_s)
    # from scipy.sparse.linalg import cg
    # S = cg(L_s, R_s)  # sparse matrix solve
    # S = np.linalg.solve(L_s,R_s)

    #get the temperature
    T = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            # Globe index
            k = j * ny + i
            # reload T
            T[i, j] = S[k]
    #calc heat flow
    hf = (T[1]-T[0])/ystp*K_layer[0]  # hf = G*K
    hf = smooth(hf,int(nx/40)*2+1)

    return x_arr,y_arr,T,hf

if __name__ == "__main__":
    # +++++++++++++++++input values++++++++++++++
    # depth_layer = [2000, 4000, 8000, 20000]  # depth of each layer in m
    layers_2D = pd.read_excel('Crust_structure_2D_test.xlsx')
    number_of_layers = int(len(layers_2D.columns)/2)
    #  from top to bottom
    layer_funcs = []
    for i in range(number_of_layers):
        x = layers_2D.iloc[:,2*i].dropna()
        y = layers_2D.iloc[:,2*i+1].dropna()
        fc_temp = interpolate.interp1d(x,y)
        layer_funcs.append(fc_temp)

    Thermal_paras = pd.read_excel('Thermal_para_2D.xlsx')

    # K_layer = [2, 3, 3, 2]  # thermal conductivity of each layer in W/mK
    # A_layer = [2.5e-6, 2e-6, 1.5e-6, 1.5e-6]  # heat production in W/m3
    K_layer = Thermal_paras['K'] # thermal conductivity of each layer in W/mK
    A_layer =  Thermal_paras['A']# heat production in W/m3
    # boundary condition
    T0 = 0  # top boundary in C
    Tb = 500  # bottom boundary in C
    # bottom_bound = Tb
    # boundary_cond = 'TEMP'
    upper_bounds = 0
    boundary_cond = 'HF'
    bottom_bd = pd.read_excel('Bottom_boundary_2D_test.xlsx')
    bottom_bounds = interpolate.interp1d(bottom_bd.iloc[:,0], bottom_bd.iloc[:, 1])
    # nz,tep = func_2D(depth_layer, k_layer, A_layer, T0, bottom_bound, boundary_cond)
    xlim = layers_2D.iloc[-1,0]  #the last number in first col
    ylim = layers_2D.iloc[:,-1][0]  #the first number in last columns

    x_arr, y_arr, T ,hf = func_2D(xlim,ylim,layer_funcs,K_layer,A_layer,upper_bounds,bottom_bounds,boundary_cond)


    fig_2d_temp = plt.figure(figsize=(10,10))

    rect1 = [0.14, 0.75, 0.77, 0.2]  # [left, bottom, width, height] the area of the figure （全部是0~1之间的数，表示比例）
    rect2 = [0.14, 0.3, 0.77, 0.4]
    rect3 = [0.14, 0.05, 0.77, 0.2]

    ax_shf = plt.axes(rect1)
    ax_shf.plot(x_arr, hf*1000)
    ax_shf.set_xlim(x_arr[0] ,x_arr[-1] )
    ax_shf.set_title('Surface Heat Flow mW/m2')

    ax_temp_field = plt.axes(rect2)
    #plot the layers
    for fc in layer_funcs:
        ax_temp_field.plot(x_arr,fc(x_arr),c='k')
    c = ax_temp_field.pcolor(x_arr, y_arr, T,cmap = 'jet')
    ax_temp_field.invert_yaxis()
    cbar_ax = fig_2d_temp.add_axes([0.92, 0.3, 0.015, 0.4])  # x,y of bottom left point, width, heigth of the colorbar
    fig_2d_temp.colorbar(c, cax=cbar_ax)
    ax_temp_field.set_title('Temperature Distribution C')
    ax_temp_field.set_ylabel('Depth km')

    ax_bhf = plt.axes(rect3)
    ax_bhf.plot(x_arr,bottom_bounds(x_arr)*1000)
    ax_bhf.set_xlim(x_arr[0],x_arr[-1])
    ax_bhf.set_title('Bottom Heat Flow mW/m2')
    ax_bhf.set_xlabel('Distance km')

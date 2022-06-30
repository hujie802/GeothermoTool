# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:54:52 2019

@author: hujie
"""


import matplotlib.pyplot as plt
from decompaction import decomp
import numpy as np
from read_file import read_file
from sedi_load import sedi_load
from paleosealevel import paleosealevel
from paleobathymetry import paleobathymetry
import pandas as pd
import os

def cal_subsidence(lay_file,phi0_sand, phi0_mud, phi0_carbon, phi0_coal, c_sand, c_mud, c_carbon, c_coal, rho_sand, rho_mud,
                       rho_carbon, rho_coal, rho_water,rho_mantle,PHI):
    #PHI is Compensation degree
    #read layer files, return dataframe
    layers = read_file(lay_file, phi0_sand, phi0_mud, phi0_carbon, phi0_coal, c_sand, c_mud, c_carbon, c_coal, rho_sand, rho_mud,
                       rho_carbon, rho_coal, rho_water)  # read file
    nlayers = decomp(layers)
    sed, t1 = sedi_load(nlayers, rho_mantle, rho_water)

    # paleo-sea-level correction
    paleos, t2 = paleosealevel(layers, rho_mantle, rho_water)

    # paleobathymetry correction
    paleob = paleobathymetry(layers)
    # cal tectonic subsudence
    Y = PHI * (sed - paleos) + paleob
    Y = np.append(Y, 0)

    subsidence = pd.DataFrame({'age':t1,'subsidence':Y})

    return subsidence


if __name__ == "__main__":
    # surface porosity
    phi0_sand = 0.6
    phi0_mud = 0.5
    phi0_carbon = 0.6
    phi0_coal = 0.9
    # coefficient  m^-1
    c_sand = 0.217e-3
    c_mud = 0.515e-3
    c_carbon = 0.22e-3
    c_coal = 0.7e-3
    # density kg/m^3
    rho_sand = 2800
    rho_mud = 2400
    rho_carbon = 2720
    rho_coal = 1800

    rho_water = 1000
    rho_mantle = 3300

    ##补偿度
    PHI = 1
    sub = cal_subsidence('sectiondata1_0.csv', phi0_sand, phi0_mud, phi0_carbon, phi0_coal, c_sand, c_mud, c_carbon, c_coal, rho_sand,
                   rho_mud,
                   rho_carbon, rho_coal, rho_water, PHI)

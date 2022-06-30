# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:27:32 2019

@author: hujie
"""

from numpy import exp
class Layer(object):
# This file is used to save the parameters of the formations.

    
    def __init__(self,formation,age_start,age_end,top,bottom,paleobathymetry,paleosealevel,
                 sand,mud,carbon,coal,
                 phi0_sand,phi0_mud,phi0_carbon,phi0_coal,
                 c_sand,c_mud,c_carbon,c_coal,
                 rho_sand,rho_mud,rho_carbon,rho_coal,rho_water):
    # formation: formation name
    # age_start: depositional starting age of this formation, Ma
    # age_end: depositional starting ending of this formation, Ma
    # top: depth of the top, m
    # bottom: depth of the bottom, m
    # sand,mud,carbon,coal: proportion of those compositions, 0-1
    # phi0_sand,phi0_mud,phi0_carbon,phi0_coal: initial porosity , 0-1
    # c_sand,c_mud,c_carbon,c_coal: Compaction coefficient, m^-1
    # rho_sand,rho_mud,rho_carbon,rho_coal,rho_water: density,kg/m3

        
        self.age_start = age_start
        self.age_end = age_end
        self.top = top
        self.bottom = bottom
        self.paleobathymetry = paleobathymetry
        self.paleosealevel = paleosealevel
        self.sand = sand
        self.formation = formation
        
        self.sand = sand
        self.mud = mud
        self.carbon = carbon
        self.coal = coal
                
        self.phi0_sand = phi0_sand
        self.phi0_mud = phi0_mud
        self.phi0_carbon = phi0_carbon
        self.phi0_coal = phi0_coal
        
        self.c_sand = c_sand
        self.c_mud = c_mud
        self.c_carbon = c_carbon
        self.c_coal = c_coal
        
        self.rho_sand = rho_sand
        self.rho_mud = rho_mud
        self.rho_carbon = rho_carbon
        self.rho_coal = rho_coal
        self.rho_water = rho_water
    
    def thickness(self):
        return self.bottom-self.top
    def phi0(self):
        return self.sand*self.phi0_sand+self.mud*self.phi0_mud+self.carbon*self.phi0_carbon+self.coal*self.phi0_coal
    def c(self):
        return self.sand*self.c_sand+self.mud*self.c_mud+self.carbon*self.c_carbon+self.coal*self.c_coal
    
    def density(self):
        phi = self.phi0()*exp(-self.c()*(self.bottom+self.top)/2)
        return (self.sand*self.rho_sand+self.mud*self.rho_mud+self.carbon*self.rho_carbon+self.coal*self.rho_coal)*(1-phi)+self.rho_water*phi
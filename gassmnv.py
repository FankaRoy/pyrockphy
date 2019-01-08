# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 05:23:16 2018

Python code by: Fanka, W. Roye T.

Ref: MATLAB code by T. Mukerji; 
     function [vp2,vs2,ro2,k2]=gassmnv(vp1,vs1,ro1,rofl1,kfl1,rofl2,kfl2,k0,\
                                          phi)
"""
#import matplotlib.pyplot as plt
import numpy as np

def Gassmnv(vp1,vs1,ro1,rofl1,kfl1,rofl2,kfl2,k0,phi):
    """
    Gassmnv(vp1,vs1,ro1,rofl1,kfl1,rofl2,kfl2,k0,phi)
    GASSMNV Gassmann's relation (velocities)
    
    Gassmann fluid substituion with velocities as input/outputs
    # =========================================================================
    #                               INPUTS
    # =========================================================================
     
     VP1, VS1, RO1: rock Vp, Vs, and density with fluid 1
     ROFL1, KFL1:   density and bulk modulus of initial fluid
     ROFL2, KFL2:   density and bulk modulus of new fluid
     K0, PHI:       mineral bulk modulus, and rock porosity
    # =========================================================================
    #                               OUTPUTS
    # =========================================================================
    Function returns a list of the following outputs:
     VP2,VS2,RO2, K2:  Vp, Vs, density, and bulk modulus of rock with new fluid
     Giving VS1=0, and mineral P-wave modulus in place of K0 does approximate 
     Gassmann calculation.
    
    *See also GASSMNK

    """

    ro2=ro1 - phi*rofl1 +phi*rofl2;
    mu1=ro1*vs1**2; k1=ro1*vp1**2-(4/3)*mu1;
    a= k1/(k0-k1) - kfl1/(phi*(k0-kfl1)) + kfl2/(phi*(k0-kfl2)); 
    k2= k0*a /(1+a); mu2=mu1;
    vp2=np.sqrt((k2+(4/3)*mu2)/ro2); 
    vs2=np.sqrt(mu2/ro2);
    
#    plt.plot(phi,vp2,'g',phi,vs2,'r')
    
    outlistgv = [vp2,vs2,ro2,k2]
    return outlistgv

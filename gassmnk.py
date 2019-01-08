# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 20:02:36 2018

@author: Fanka, W. Roye T.

Ref: MATLAB code by T. Mukerji; 
     function k2=gassmnk(k1,kfl1,kfl2,k0,phi)
"""
import matplotlib.pyplot as plt

def Gassmnk(k1,kfl1,kfl2,k0,phi):
    """
    Gassmnk(k1,kfl1,kfl2,k0,phi)
    GASSMNK Gassmann's relation (moduli)
    
    Fluid substitution using Gassmann equation
     k1:  original bulk modulus of rock saturated with fluid of bulk modulus 
          kfl1
     k2:  rock bulk modulus with new fluid of bulk modulus KFL2
     k0:  mineral modulus
     phi: porosity
    
     Giving the P-wave modulus in place of the bulk modulus does the 
     approximate Gassmann calculation.
    
    *See also GASSMNV
    """

    a= k1/(k0-k1) - kfl1/(phi*(k0-kfl1)) + kfl2/(phi*(k0-kfl2)); 
    k2= (k0*a /(1+a))*(phi!=0) + k1*(phi==0); 
    
    plt.plot(phi,k2,'b')
    
    return k2

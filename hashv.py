# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 18:25:08 2018

@author: Fanka, W. Roye T.

Ref: MATLAB code by T. Mukerji; 
    function [vpu,vpl,vsu,vsl,por]=hashv(vp1,vs1,ro1,vp2,vs2,ro2)
"""
import matplotlib.pyplot as plt
import numpy as np

outlistv = [] # initial 'global' list of outputs.
def Hashv(vp1,vs1,ro1,vp2,vs2,ro2):
    """
    Function computes velocities for the Hashin-Shtrikman 
    upper and lower bound curves.
    Without input arguments, plots default velocity bounds as a function of 
    fraction of material 2.
    # =========================================================================
    #                           INPUTS
    # =========================================================================
    Takes the following inputs:
        vp1,vs1,vp2,vs2:   Velocities of the two constituents
        ro1,ro2:           Density of the two constituents
    # =========================================================================
    #                           OUTPUTS
    # =========================================================================
    Function returns a list of outputs:
        [VPU,VPL,VSU,VSL,POR,RHO]
        VPU,VPL,VSU,VSL:       Upper and lower bounds on velocities
        POR:                   Volume fraction of material 2
        RHO:                Bulk density
    # =========================================================================
    #                           ASSUMPTION  
    # =========================================================================
    Assumes material 1 has higher velocity than material 2. If not, then upper 
    and lower bounds should be interchanged in the output.
    
    *See also HASH
    """
    mu1=ro1*vs1**2; mu2=ro2*vs2**2;
    k1=ro1*vp1**2-(4/3)*mu1; k2=ro2*vp2**2-(4/3)*mu2;
    por=np.linspace(0,1,100); por[0]=1e-7;
    ku=k2+(1.-por)*(k1-k2)*(k2+4.*mu1/3.)/(k2+4.*(mu1/3.)+por*(k1-k2));
    kl=k2+(1.-por)*(k1-k2)*(k2+4.*mu2/3.)/(k2+4.*(mu2/3.)+por*(k1-k2));
    fgu=mu1*(9.*k1+8.*mu1)/(6.*(k1+2.*mu1));
    fgl=mu2*(9.*k2+8.*mu2)/(6.*(k2+2.*mu2));
    gu=mu2+(mu1-mu2)*(1.-por)*(mu2+fgu)/(mu2+fgu+por*(mu1-mu2));
    gl=mu2+(mu1-mu2)*(1.-por)*(mu2+fgl)/(mu2+fgl+por*(mu1-mu2));

#    gl=mu1*(1.-por)/(1.+por*mu1);
    ro=(1.-por)*ro1+por*ro2; rho=ro;
    vpu=np.sqrt((ku+(4.*gu/3.))/ro);
    vpl=np.sqrt((kl+(4.*gl/3.))/ro);
    vsu=np.sqrt(gu/ro);
    vsl=np.sqrt(gl/ro);

    plt.plot(por,vpu,'-g',por,vpl,'-g',por,vsu,'--c',por,vsl,'--c')
    
    outlistv = [vpu,vpl,vsu,vsl,por,rho]
    
    return outlistv


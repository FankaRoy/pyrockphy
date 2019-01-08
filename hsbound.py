# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 20:25:33 2018

@author: Fanka
"""
import matplotlib.pyplot as plt
import numpy as np

def HSbound(k1,mu1,k2,mu2,por):
    global KnG
    """
    HSBOUND Hashin-Shtrikman bound 
    example
    por=np.linspace(0,.4,100); # quartz=1, shale=2
    HSbound(38,44,21.7,6.67,por)
    """

    ku=k2+(1.-por)*(k1-k2)*(k2+4.*mu1/3.)/(k2+4.*(mu1/3.)+por*(k1-k2));
    kl=k2+(1.-por)*(k1-k2)*(k2+4.*mu2/3.)/(k2+4.*(mu2/3.)+por*(k1-k2));
    fgu=mu1*(9.*k1+8.*mu1)/(6.*(k1+2.*mu1));
    fgl=mu2*(9.*k2+8.*mu2)/(6.*(k2+2.*mu2));
    gu=mu2+(mu1-mu2)*(1.-por)*(mu2+fgu)/(mu2+fgu+por*(mu1-mu2));
    gl=mu2+(mu1-mu2)*(1.-por)*(mu2+fgl)/(mu2+fgl+por*(mu1-mu2));
    
    plt.plot(por,ku,'k',por,kl,'b',por,gu,'r',por,gl,'m')


    
#    ku=list(ku); kl=list(kl); gu=list(gu); gl=list(gl); 
    KnG = [ku, kl, gu, gl]
    return np.array(KnG)

#por=np.linspace(0,1,1000);
#HSbound(38,44,21.7,6.67,por)

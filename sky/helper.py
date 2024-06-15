# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 03:02:04 2016

@author: shimura
"""

def raman_shift(x_scattered, x_raman_pump):
        """
        Parameters
        ----------
        x_scattered  :  散乱波長 [nm]
        x_Rpump      :  Raman pump波長 [nm]
        
        Return
        ------
        q  :  Raman shift wavenumber [cm-1]
        if q < 0: Stokes, if q > 0: anti-Stokes scattering
        """
        q = ( 1./x_scattered - 1./x_raman_pump )*1.0e7
        return q
     
def inverse_raman_shift(q, x_raman):
    """
    raman_shift の逆変換（q --> x_scattered）
    Rdatで利用する。
    """
    x_scattered = 1./( q*1.0e-7 + 1./x_raman)
    return x_scattered
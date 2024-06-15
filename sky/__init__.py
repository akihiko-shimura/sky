# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 02:56:37 2016

@author: A.Shimura

This module is useful functions for treating sky-format data in Python.
Sky format is defined in the femtosecond absorption measurement program, 
fs_kHz.exe (made by Masayuki Yoshizawa).
Python3 version updated in 2020/08/17

"""
__version__ = "0.1.0"

from . import classes, helper
from .classes import *

__all__ = ["classes", "helper"]

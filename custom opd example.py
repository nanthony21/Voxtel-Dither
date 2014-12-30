# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 17:51:37 2014

@author: Nick
"""

import numpy as np
import matplotlib.pyplot as plt
import Dither

a=np.zeros((335,335))
r=1*a
for i in range(335):
    for j in range(335):
        r[i,j]=2.1215*np.sqrt((i-335/2)**2+(j-335/2)**2)/(335/2)
r[r>2.12]=None

n=1.5438-.0083116314*r**2


d=Dither.Dither(0,15,6,1.5064,1.5438,idealn=n)
d.calc3d()
d.savelayers('New folder')
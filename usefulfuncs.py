# -*- coding: utf-8 -*-
"""
Created on Fri Aug 08 14:11:51 2014

@author: CHANGME
"""
from __future__ import division

def optimalheight(opd,vheight,n1,n2):
    '''
    Given the desired optical path difference (opd), the height of each voxel (vheight),
    and the two indices of refraction, this function with tell you the minimum needed layers
    that are needed to acheive the desired opd.
    '''
    if n1>n2:
        raise ValueError('n2 must be greater than n1')
    optheight=opd/(vheight*(n2-n1))
    return optheight
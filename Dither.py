# -*- coding: utf-8 -*-
"""
Created on Fri Aug 01 15:43:31 2014

@author: Nick
"""

from __future__ import division
import numpy as np, matplotlib.pyplot as plt
from scipy import misc
import bisect



class Dither:
    
    def __init__(self,ideal,height,pheight,n1,n2):
        self.height=height
        
        self.ideal=ideal-np.nanmin(ideal)
        self.ypix=ideal.shape[0]
        self.xpix=ideal.shape[1]
        self.pheight=pheight
        self.n1=n1
        self.n2=n2
        if n1>n2:
            raise ValueError('N2 must be greater than N1')
        self.possiblevalues=[]
        for i in xrange(self.height+1):
            self.possiblevalues.append((i*self.n2+(self.height-i)*self.n1)/(self.height))
        self.idealn=self.ideal/(self.height*self.pheight)+self.n1
        self.__floyddither()
    
        
     

                   
    def __floyddither(self):
        '''
        This takes the ideal OPD and quantizes it and dithers it based on the number available
        opds based on n1, n2, and the number of layers. It uses the Floyd Steinberg algorithm.
        Based on psuedo code found at http://en.wikipedia.org/wiki/Floyd%E2%80%93Steinberg_dithering
        '''
        self.dithered=1*self.idealn
        
        mids = [(self.possiblevalues[i] + self.possiblevalues[i + 1]) / 2.0 for i in xrange(len(self.possiblevalues) - 1)]       
        
        for i in xrange(self.ypix):
            for j in xrange(self.xpix):
                if np.isnan(self.idealn[i,j]):
                    self.dithered[i,j]=None
                else:
                    ind = bisect.bisect_right(mids, self.dithered[i,j])
                    new=self.possiblevalues[ind]                    
                    qerror=self.dithered[i,j]-new
                    self.dithered[i,j]=1*new
                    try:
                        self.dithered[i,j+1]=self.dithered[i,j+1]+qerror*0.4375
                    except IndexError:
                        pass                   
                    try:
                        self.dithered[i+1,j-1]=self.dithered[i+1,j-1]+qerror*0.1875
                    except IndexError:
                        pass                    
                    try:
                        self.dithered[i+1,j]=self.dithered[i+1,j]+qerror*0.3125
                    except IndexError:
                        pass                    
                    try:
                        self.dithered[i+1,j+1]=self.dithered[i+1,j+1]+qerror*0.0625
                    except IndexError:
                        pass
    
    def calc3d(self):         
              

        self.x1=np.uint8(np.round(self.height*(self.dithered-self.n2)/(self.n1-self.n2)))
    
        self.data=np.uint8(self.x1[...,np.newaxis]>np.arange(self.height))       
        
        self.data=self.data.swapaxes(-1,-1)
        shp=self.data.shape[:-1]
        for ndx in np.ndindex(shp):
            np.random.shuffle(self.data[ndx])
    
        
    
    def printlayers(self):
        #shows each layers image.
        for i in range(self.height):
            plt.figure()    
            plt.imshow(self.data[:,:,i],interpolation='none',cmap='gray')
        plt.show()
    
    def savelayers(self, directory):
        n1Data = self.data
        n2Data = self.data==False        
                        
        for i in range(self.height):
            misc.imsave(directory + '\layer %d (n1).bmp'%(i+1),n1Data[:,:,i].astype(int))
            misc.imsave(directory + '\layer %d (n2).bmp'%(i+1),n2Data[:,:,i].astype(int))
     
    def image(self):
         '''
         plots an images showing the original OPD, the quantized OPD, and the dithered OPD.
         '''

         fig=plt.figure()
         ax=fig.add_subplot(2,2,1)
         ax3=fig.add_subplot(2,2,3)
         ax.imshow(self.idealn)
         ax3.imshow(self.dithered)
         

        
        
        
def examplegauss(res,sigma):
    gauss=np.zeros((res,res))
    for i in range(res):
        for j in range(res):
            r2=(i-res/2)**2+(j-res/2)**2
            gauss[i,j]=np.exp(-r2/(2*sigma**2))
    return gauss
        
      


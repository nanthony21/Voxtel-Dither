# -*- coding: utf-8 -*-
"""
Created on Fri Aug 01 15:43:31 2014

@author: Nick
"""

from __future__ import division
import scipy as sp, numpy as np, matplotlib.pyplot as plt
from scipy import misc
from random import randrange
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
    
    def calc3d(self,lim=None):
        '''
        This function performs the more time consuming 3d dimensional calculations.
        Since they may not always be needed they are not automatically performed.
        This function will need to be run before printlayers() or savelayer() can be
        performed or the "data" array can be accesed.
        '''
        self.lim=lim
        self.__possibilities(self.dithered)
        self.__selection1()
        
    def calc3dfast(self):
        '''
        This function determines x1, the number of n1 dots required for each pixel.
        It then creates the 3d array, data, that contains either 1 or 0 where 1 represents an n1 voxel
        and a 0 represents a n2 voxel. All the n1 voxels are stacked at the bottom with all the n2 voxels
        stacked on top.
        '''
        self.x1=np.uint8(np.round(self.height*(self.dithered-self.n2)/(self.n1-self.n2)))
        self.data=sp.zeros([self.ypix,self.xpix,self.height],dtype=np.uint8)
        for i in xrange(self.height):
            self.data[:,:,i][self.x1>=i] = 1

     
    def __quantize(self):
        '''       
        given "ideal", the numpy array representing the ideal image in the transverse plane,
        "height", the number of printed layers (bit depth), and "n1" and "n2" the two
        possible values for each voxel, this function creates a quantized version of ideal.
        '''        
        self.quantized=1*self.idealn
        for i in xrange(self.ypix):
            for j in xrange(self.xpix):
                if np.isnan(self.idealn[i,j]):
                    self.quantized[i,j]=None
                else:
                    self.quantized[i,j]=self.__quantizer(self.idealn[i,j],self.possiblevalues)
    def __floyddither(self):
        '''
        This takes the ideal OPD and quantizes it and dithers it based on the number available
        opds based on n1, n2, and the number of layers. It uses the Floyd Steinberg algorithm.
        Based on psuedo code found at http://en.wikipedia.org/wiki/Floyd%E2%80%93Steinberg_dithering
        '''
        self.dithered=1*self.idealn
        for i in xrange(self.ypix):
            for j in xrange(self.xpix):
                if np.isnan(self.idealn[i,j]):
                    self.dithered[i,j]=None
                else:
                    old=1*self.dithered[i,j]            
                    new=self.__quantizer(self.dithered[i,j],self.possiblevalues)
                    self.dithered[i,j]=1*new            
                    qerror=old-new
                    if j<=(self.xpix-2):
                        self.dithered[i,j+1]=self.dithered[i,j+1]+qerror*7/16
                    if i<=(self.ypix-2) and j>=(1):
                        self.dithered[i+1,j-1]=self.dithered[i+1,j-1]+qerror*3/16
                    if i<=(self.ypix-2):
                        self.dithered[i+1,j]=self.dithered[i+1,j]+qerror*5/16
                    if i<=(self.ypix-2) and j<=(self.xpix-2):
                        self.dithered[i+1,j+1]=self.dithered[i+1,j+1]+qerror*1/16
                    
    
    def __possibilities(self,inputarray):
        def combos():
            '''
            Given one argument, height, which is essentially bitdepth, this function returns two things.
            It returns "combos", a list of lists of list combinations.
            It also returns "numadj", the number of adjacent identical bits in the corresponding combo.
            
            These are both automatically sorted so that the combos with the lowest numadj come first.
            They both follow the same indexing scheme. If the first index is "i", then you have selecte
            the list of combos for which there are "i" HIGH bits. The next index selects the particular
            combo that you are interested in. Presumably you will prefer the combos with the lower second
            indexes because they have less adjacent like bits and are therefore better
            as far as dithering is concerned.
            
            The optional "lim" argument throws away unwanted combinations. For example, if "lim"=2 then only the combinations
            with the two lowest numadj in the set are returned. 
            
            "lim" must be less than "height" and they must both be integers.
            '''
        
            def comboadjnum(combo):
                n=[]
                for i in range(len(combo)):  
                    n.append([])
                    for j in range(len(combo[i])):        
                        b=2  
                        a=0
                        for k in combo[i][j]:
                            if k==b:
                                a+=1
                            b=k
                        n[i].append(a)
                return n
            
            def combo():
                def binarize(num,length):
                    a=bin(num)[2:]
                    b=len(a)
                    c=sp.zeros(length,dtype=np.uint8)
                    for i in range(b):
                        c[-(i+1)]=a[-(i+1)]
                    return list(c)
                dotnum=[[]for x in xrange(self.height+1)]
                for j in xrange(2**self.height):
                    bina=binarize(j,self.height)
                    dotnum[sum(bina)].append(bina)
                return dotnum
            
            
            #this sorts the combos in order of numadj
            combos=[]   
            numadj=[]
            c=combo()
            n=comboadjnum(c)
            for i in range(len(c)):    
                b=sorted(zip(n[i],c[i]))
                combos.append([point[1] for point in b])
                numadj.append([point[0] for point in b])
            
            if self.lim:    
                combolim=[[] for x in xrange(len(combos))]
                numadjlim=[[] for x in xrange(len(combos))]
                for i in range(len(combos)):
                    limnum=(list(set(numadj[i]))+[self.height-1]*(self.height-self.lim))[self.lim-1]
                
                    for j in xrange(len(combos[i])):
                         if numadj[i][j]<=limnum:
                             combolim[i].append(combos[i][j])
                             numadjlim[i].append(numadj[i][j])
                combos=combolim
                numadj=numadjlim
            return combos,numadj        
        
        
        
        
        #determines x1, the number of n1 dots required for each pixel.
        self.x1=np.uint8(np.round(self.height*(inputarray-self.n2)/(self.n1-self.n2)))
        co=combos()[0]
    
        #generates 2d array of all possible combo lists
        self.possibilities=sp.zeros(inputarray.shape,dtype='O')
        for y in xrange(self.ypix): 
            for x in xrange(self.xpix):
                self.possibilities[y,x]=co[self.x1[y,x]]
    
    def __selection1(self):
        '''
        This creates  the 3d array, Data, by selecting out one of the 1d possibility
        arrays to be used for each of the values of the 2d array, Dithered.
        '''
        self.data=sp.zeros([self.ypix,self.xpix,self.height],dtype=np.uint8)
        for y in xrange(self.ypix): 
            for x in xrange(self.xpix):
                self.data[y,x]=self.possibilities[y,x][randrange(len(self.possibilities[y,x]))]
    
    
    def printlayers(self):
        #shows each layers image.
        for i in range(self.height):
            plt.figure()    
            plt.imshow(self.data[:,:,i],interpolation='none',cmap='gray')
            plt.show()
    
    def savelayers(self, directory):
        n1Data = self.data
        n2Data = self.data==False        
        l = len(n1Data)        
                        
        for i in range(self.height):
            misc.imsave(directory + '\layer %d (n1).bmp'%(i+1),n1Data[:,:,i].astype(int))
            misc.imsave(directory + '\layer %d (n2).bmp'%(i+1),n2Data[:,:,i].astype(int))
     
    def image(self):
         '''
         plots an images showing the original OPD, the quantized OPD, and the dithered OPD.
         '''
         self.__quantize()
         fig=plt.figure()
         ax=fig.add_subplot(2,2,1)
         ax2=fig.add_subplot(2,2,2)
         ax3=fig.add_subplot(2,2,3)
         ax.imshow(self.idealn)
         ax2.imshow(self.quantized)
         ax3.imshow(self.dithered)
        
    def qerrorimage(self):
        '''
        plots an image of the the original OPD minus the dithered OPD, it is a
        representation of the quantization error.
        '''
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        p=ax.imshow(self.idealn-self.dithered,interpolation='none')
        plt.colorbar(p,ax=ax)
        
    def __quantizer(self,num, quant):
        mids = [(quant[i] + quant[i + 1]) / 2.0
                for i in xrange(len(quant) - 1)]
        ind = bisect.bisect_right(mids, num)
    
        return quant[ind]
      


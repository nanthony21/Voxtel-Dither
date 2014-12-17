from __future__ import division
import numpy as np, scipy as sp
import matplotlib.pyplot as plt

class Zernike:
    '''
    Given a list, "cofs", of the zernike coefficients the wavelength and the
    desired resolution in the X and Y direction this class creates a variety of useful outputs.
    
    For a given instance of this class, z, you have access to:
    z.x,z.y,and z.r: the relevant coordinate arrays.
    z.wf: the wavefront in units of wavelengths and centered at 0.
    z.opd: the optical path difference needed to create the desired wavefront.
    This is in the same units that were input for wavelength.
    
    z.plot3d() will plot z.opd in 3d
    
    z.image() will create a colored image representing z.opd
    '''

    def delta(a,b):
        if a == b:
            return 1
        else:
            return 0
    
    def __init__(self,cofs,wavelength,res,normalizationratio,correcting=False,fringe=False):
        self.xres = self.yres = res                
        X = np.linspace(-1,1,self.xres)
        Y = np.linspace(-1,1,self.yres)
        self.cofs = cofs
        if len(self.cofs) < 37:
            self.cofs = np.append(self.cofs,np.zeros(37-len(self.cofs)))
        self.wavelength = wavelength
        self.x, self.y =np.meshgrid(X,Y)
        self.r =np.sqrt(self.x**2+self.y**2)
        self.theta = np.arctan2(self.y,self.x)
        self.wf=sp.zeros((self.yres,self.xres))
        if normalizationratio<1:
            raise ValueError('normalizedradius may not be less than printedradius')
        for i in range(37):
            if self.cofs[i]!=0 and fringe==False:
               self.wf+= self.cofs[i]*getattr(self,'z%g'%i)(self.r/normalizationratio,self.theta)
            if self.cofs[i]!=0 and fringe==True:
               self.wf+= self.cofs[i]*getattr(self,'zf%g'%i)(self.r/normalizationratio,self.theta)
               
        self.wf[self.r > 1] = None
        
        if correcting==True:
            self.wf=-1*self.wf
        
        nanmin=np.nanmin(self.wf)
        self.opd=((self.wf-nanmin)*self.wavelength)
        self.params = [self.opd, self.x, self.y]
          
    
    
    def z0(self,r,theta):
        return 1
    def z1(self,r,theta):
        return np.sqrt (4)*r*np.cos(theta)
    def z2(self,r,theta):
        return np.sqrt (4)*r*np.sin(theta)
    def z3(self,r,theta):
        return np.sqrt (3)*(2*r**2-1)
    def z4(self,r,theta):
        return np.sqrt (6)*r**2*np.sin(2*theta)
    def z5(self,r,theta):
        return np.sqrt (6)*r**2*np.cos(2*theta)
    def z6(self,r,theta):
        return np.sqrt (8)*(3*r**3-2*r)*np.sin(theta)
    def z7(self,r,theta):
        return np.sqrt (8)*(3*r**2-2*r)*np.cos(theta)
    def z8(self,r,theta):
        return np.sqrt (8)*r**3*np.sin(3*theta)
    def z9(self,r,theta):
        return np.sqrt (8)*r**3*np.cos(3*theta)
    def z10(self,r,theta):
        return np.sqrt (5)*(6*r**4-6*r**2+1)
    def z11(self,r,theta):
        return np.sqrt (10)*(4*r**4-3*r**2)*np.cos(2*theta)
    def z12(self,r,theta):    
        return np.sqrt (10)*(4*r**4-3*r**2)*np.sin(2*theta)
    def z13(self,r,theta):
        return np.sqrt (10)*r**4*np.cos(4*theta)
    def z14(self,r,theta):
        return np.sqrt (10)*r**4*np.sin(4*theta)
    def z15(self,r,theta):
        return np.sqrt (12)*(10*r**5-12*r**3+3*r)*np.cos(theta)
    def z16(self,r,theta):
        return np.sqrt (12)*(10*r**5-12*r**3+3*r)*np.cos(theta)
    def z17(self,r,theta):
        return np.sqrt (12)*(5*r**5-4*r**3)*np.cos(3*theta)
    def z18(self,r,theta):
        return np.sqrt (12)*(5*r**5-4*r**3)*np.sin(3*theta)
    def z19(self,r,theta):
        return np.sqrt (12)*(r**5)*np.cos(5*theta)
    def z20(self,r,theta):
        return np.sqrt (12)*(r**5)*np.sin(5*theta)
    def z21(self,r,theta):
        return np.sqrt (7)*(20*r**6-30*r**4+12*r**2-1)
    def z22(self,r,theta):
        return np.sqrt (14)*(15*r**6-20*r**4+6*r**2)*np.sin(2*theta)
    def z23(self,r,theta):
        return np.sqrt (14)*(15*r**6-20*r**4+6*r**2)*np.cos(2*theta)
    def z24(self,r,theta):
        return np.sqrt (14)*(6*r**6-5*r**4)*np.sin(4*theta)
    def z25(self,r,theta):
        return np.sqrt (14)*(6*r**6-5*r**4)*np.cos(4*theta)
    def z26(self,r,theta):
        return np.sqrt (14)*(r**6)*np.sin(6*theta)
    def z27(self,r,theta):
        return np.sqrt (14)*(r**6)*np.cos(6*theta)
    def z28(self,r,theta):
        return np.sqrt (16)*(35*r**7-60*r**5+30*r**3-4*r)*np.sin(theta)
    def z29(self,r,theta):
        return np.sqrt (16)*(35*r**7-60*r**5+30*r**3-4*r)*np.cos(theta)
    def z30(self,r,theta):
        return np.sqrt (16)*(21*r**7-30*r**5+10*r**3)*np.sin(3*theta)
    def z31(self,r,theta):
        return np.sqrt (16)*(21*r**7-30*r**5+10*r**3)*np.cos(3*theta)
    def z32(self,r,theta):
        return np.sqrt (16)*(7*r**7-6*r**5)*np.sin(5*theta)
    def z33(self,r,theta):
        return np.sqrt (16)*(7*r**7-6*r**5)*np.cos(5*theta)
    def z34(self,r,theta):
        return np.sqrt (16)*(r**7)*np.sin(7*theta)
    def z35(self,r,theta):
        return np.sqrt (16)*(r**7)*np.cos(7*theta)
    def z36(self,r,theta):
        return np.sqrt (9)*(70*r**8-140*r**6+90*r**4-20*r**2+1)
        
    def zf0(self,r,theta):
        return 1    
    def zf1(self,r,theta):
        return r*np.cos(theta)
    def zf2(self,r,theta):
        return r*np.sin(theta)        
    def zf3(self,r,theta):
        return 2*r**2-1
    def zf4(self,r,theta):
        return r**2*np.cos(2*theta)
    def zf5(self,r,theta):
        return r**2*np.sin(2*theta)  
    def zf6(self,r,theta):
        return (3*r**2-2)*r*np.cos(theta)
    def zf7(self,r,theta):
        return (3*r**2-2)*r*np.sin(theta)
    def zf8(self,r,theta):
        return 6*r**4-6*r**2+1
    def zf9(self,r,theta):
        return r**3*np.cos(3*theta)
    def zf10(self,r,theta):
        return r**3*np.sin(3*theta)
    def zf11(self,r,theta):
        return (4*r**2-3)*r**2*np.cos(2*theta)
    def zf12(self,r,theta):
        return (4*r**2-3)*r**2*np.sin(2*theta)
    def zf13(self,r,theta):
        return (10*r**4-12*r**2+3)*r*np.cos(theta)
    def zf14(self,r,theta):
        return (10*r**4-12*r**2+3)*r*np.sin(theta)
    def zf15(self,r,theta):
        return (20*r**6-30*r**4+12*r**2-1)
    def zf16(self,r,theta):
        return r**4*np.cos(4*theta)
    def zf17(self,r,theta):
        return r**4*np.sin(4*theta)
    def zf18(self,r,theta):
        return (5*r**2-4)*r**3*np.cos(3*theta)
    def zf19(self,r,theta):
        return (5*r**2-4)*r**3*np.sin(3*theta)
    def zf20(self,r,theta):
        return (15*r**4-20*r**2+6)*r**2*np.cos(2*theta)
    def zf21(self,r,theta):
        return (15*r**4-20*r**2+6)*r**2*np.sin(2*theta)
    def zf22(self,r,theta):
        return (35*r**6-60*r**4+30*r**2-4)*r*np.cos(theta)
    def zf23(self,r,theta):
        return (35*r**6-60*r**4+30*r**2-4)*r*np.sin(theta)
    def zf24(self,r,theta):
        return 70*r**8-140*r**6+90*r**4-20*r**2+1
    
        
        
    def image(self):
        plt.figure()
        plt.imshow(self.opd,cmap='jet',extent=[self.x.min(), self.x.max(), self.y.min(), self.y.max()])
        plt.colorbar()
        plt.show()
             

# -*- coding: utf-8 -*-
"""
Created on 21 Nov 2017

@author: thanagno

script for MFD & barycentric position of the slit 

"""

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from scipy import optimize
from scipy import ndimage


class gen_analysis(object):
    '''
    Class that calculates the Gaussian fit for each inputed image

    Params: 
        file (np.array): image file of 3D datacube

    '''
    def __init__(self):
        self.file = []
#
    def gauss(self, x, a, x0, sigma):
        ''' calculates a 1d gaussian, given an input array (x) and constants a,x0,sigma
            a = height of max
            x0 = centre
            sigma = width (standard deviation)'''
        a, x0, sigma = float(a), float(x0), float(sigma)
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
#
    def gauss_fit(self,file):#, start, end):
    #------------------read fits--------------------------
        a = file#[:,start:end]
    #    a = slit[:, 80:260]
        row = a.shape[1]
        x = np.arange(a.shape[0])
        wm = np.zeros(a.shape[1])
        com = np.zeros(a.shape[1])
        std = np.zeros(a.shape[1])
        for i in np.arange(row):
            y = a[:,i].copy()
            #------------------fit--------------------------
            popt, pcov = optimize.curve_fit(self.gauss, x, y, p0=[y.max(), y.shape[0]/2, 1], maxfev=100000)
            #------------------Generating a Gauss from the fit parameters--------------------------
            x1 = np.linspace(x.min(), x.max(), 1000)
            fitted_y = self.gauss(x1, popt[0], popt[1], popt[2])
            # plt.plot(x,fitted_y)
            # plt.plot(x, y, 'o')
            # plt.show()
            wm[i] = np.where(fitted_y>=fitted_y.max()*1./np.exp(2))[0].shape[0]
            com[i] = popt[1]
    #        # std[i] = np.sqrt(np.diag(pcov).mean())
    #        std[i] = np.std(fitted_y)
        return wm, com
#
    def identify_objects(self,image,threshold):
        ''' taken from the internet        
        Returns a mask of numbers, with the lowest value being the background (0), then 1 and 2 etc...
        http://stackoverflow.com/questions/5298884/finding-number-of-colored-shapes-from-picture-using-python
        '''    
    # smooth the image (to remove small objects); set the threshold
        global number
        number = 0
        dnaf = ndimage.gaussian_filter(image, 2)
        T  = threshold    #MAKE SAME AS DAVE!
        # find connected components
        labeled, nr_objects = ndimage.label(dnaf > T) # `dna[:,:,0]>T` for red-dot case
        #print "Number of objects is %d " % nr_objects
        #REMOVE SMALL/WRONG ONES....
        labeled = np.where(labeled > 0,1,0)
        labeled2 = labeled.copy()
        return labeled2

#------soapy------
path = '/home/tanagnos/Desktop/Ph.D/sim_dicer/slit_analysis/'
file = fits.getdata(path+'slit_400ms.fits')
file = np.transpose(file,(0,2,1))
file = file[:,88:-88,140:-155]

wm = np.zeros((file.shape[0], file.shape[2]))
com = np.zeros((file.shape[0], file.shape[2]))

for i in range(file.shape[0]):
    print i
    a = file[i].copy()
    f_an = gen_analysis()
    labeled = f_an.identify_objects(a, np.mean(a))
    slit = np.where(labeled == 1, a, 0)
#   check if slit file makes sense!!!!
    wm[i], com[i] = f_an.gauss_fit(slit)


#-----[see above gauss_fit x1 = np.linspace(x.min(), x.max(), 1000)]---101 pxl spaning 50.5 μm, 1000-->0.0505
wm1 = wm*0.0505
com1 = (com-com.mean())*0.0505

#------------------------------
com = np.zeros((file.shape[0], file.shape[1]))
cm = np.zeros(file.shape[0])
for i in range(file.shape[0]):
    print i
    a = file[i].copy()
    f_an = gen_analysis()
    labeled = f_an.identify_objects(a, np.mean(a))
    fiber = np.where(labeled == 1, a, 0)
    com[i] = np.sum(fiber, axis = 1)
    cm[i] = ndimage.measurements.center_of_mass(com[i])[0]

c = cm*0.5
var = (c.max() - c.min())/2
var = 100*var/8.4

a = c-(c.max()+c.min())/2
a = a*8.4/1000.
vara = (a.max() - a.min())/2    # vara = 0.00019613519274906679
#----------------------------------------
import matplotlib
font = {'family' : 'normal',
   #'weight' : 'bold',
       'size'   : 13}

matplotlib.rc('font', **font)
import matplotlib
matplotlib.rcParams.update({'font.family': 'sans-serif'})

fig = plt.figure(0)
ax1=fig.add_subplot(3,1,1)
ax2=fig.add_subplot(3,1,2)
ax3=fig.add_subplot(3,1,3)

pos1 = ax1.imshow(np.mean(file,axis=0)/np.mean(file,axis=0).max(), cmap='jet', vmin=0, vmax=1)
ax2.plot(np.mean(wm1, axis=0), '-', color='blue', ms=4)
ax2.errorbar(np.arange(582),np.mean(wm1, axis=0), yerr=np.std(wm1,axis=0), ecolor='r', capsize=2, errorevery=10, elinewidth=1)
ax3.plot(a, 'o', color='blue', ms=4)
ax2.set_xlim([0,291])

ax1.set_ylabel('Position across \nthe slit ($\mathrm{\mu}$m)')
ax1.set_xlabel('Position along the slit ($\mathrm{\mu}$m)')
ax2.set_xlabel('Position along the slit ($\mathrm{\mu}$m)')
ax2.set_ylabel('Mode Field \nDiameter \n($1/e^{2}, \mathrm{\mu}$m)')
ax3.set_ylabel('Barycentric shift \nacross the slit\n ($d$/1000)')
ax3.set_xlabel('Frame')

ax1.set_yticklabels([3, 25,-25])
ax1.set_xticklabels([5, 0, 50, 100, 150, 200, 250])
ax2.set_xticklabels([0, 50, 100, 150, 200, 250])

cbaxes1 = fig.add_axes([0.92, 0.76, 0.015, 0.2])
cb1 = fig.colorbar(pos1, cax=cbaxes1, ax=ax1, orientation="vertical")

plt.subplots_adjust(top=0.95,bottom=0.11,left=0.25,right=0.915,hspace=0.77,wspace=0.2)
# plt.savefig('slit.pdf', format='pdf')

plt.show()


# a = slit[88:-88,140:-155]
# row = a.shape[1]
# x = np.arange(a.shape[0])
# wm = np.zeros(a.shape[1])
# com = np.zeros(a.shape[1])
# std = np.zeros(a.shape[1])
# #for i in np.arange(row):
# y = a[:,100].copy()
# popt, pcov = optimize.curve_fit(gauss, x, y, p0=(y.max(), a.shape[0]/2, 2))
# x1 = np.linspace(x.min(), x.max(), 1000)
# fitted_y = gauss(x1, popt[0], popt[1], popt[2])
# plt.plot(x1,fitted_y, 'x')
# plt.plot(x, y, 'o')
# plt.show()
# w = np.where(fitted_y>=fitted_y.max()*1./np.exp(2))[0].shape[0]
# com = popt[1]
# sigma = [pcov[0,0], pcov[1,1], pcov[2,2] ]


# fig = plt.figure(0)
# ax1=fig.add_subplot(3,1,1)
# ax2=fig.add_subplot(3,1,2)
# ax3=fig.add_subplot(3,1,3)

# pos1 = ax1.imshow(np.mean(file,axis=0)/np.mean(file,axis=0).max(), cmap='jet', vmin=0, vmax=1)
# ax2.plot(np.mean(wm1, axis=0), '-', color='blue', ms=4)
# ax3.plot(np.mean(com1,axis=0), '-', color='blue')
# ax2.errorbar(np.arange(582),np.mean(wm1, axis=0), yerr=np.std(wm1,axis=0), ecolor='r', capsize=2, errorevery=10, elinewidth=1)
# ax3.errorbar(np.arange(582),np.mean(com1, axis=0), yerr=np.std(com1,axis=0), ecolor='r', capsize=2, errorevery=10, elinewidth=1)

# ax2.set_xlim([0,291])
# ax3.set_xlim([0,291])

# ax1.set_ylabel('Position across \nthe slit ($\mathrm{\mu}$m)')
# ax2.set_ylabel('Mode Field \nDiameter ($1/e^{2}, \mathrm{\mu}$m)')
# ax3.set_ylabel('Barycentre \n Position ($\mathrm{\mu}$m)')
# ax3.set_xlabel('Position along the slit ($\mathrm{\mu}$m)')
# ax1.set_yticklabels([3, 25,-25])
# ax1.set_xticklabels([5, 0, 50, 100, 150, 200, 250])
# ax2.set_xticklabels([0, 50, 100, 150, 200, 250])
# ax3.set_xticklabels([0, 50, 100, 150, 200, 250])

# cbaxes1 = fig.add_axes([0.92, 0.65, 0.015, 0.25])
# cb1 = fig.colorbar(pos1, cax=cbaxes1, ax=ax1, orientation="vertical")


# plt.subplots_adjust(top=0.88, bottom=0.11, left=0.17, right=0.915, hspace=0.2, wspace=0.2)
# #plt.savefig('slit.eps', format='eps')

# plt.show()


# SLIT = fits.getdata("cl.fits")[:,150:190,:]
# a = np.mean(SLIT, axis=0)

# #----------identify obj by sigma clipping-------
# labeled = identify_objects(a, np.mean(a))
# slit = np.where(labeled == 1, a, 0)

# plt.imshow(slit[32:49, 70:275])
# plt.show()

# wm, com = gauss_fit(slit, 84, 260)

# #----Canary data-------
# mag = 17.74674 #magnification of slit in camera
# pxl = 30 #pixel size in μm
# wm = wm*pxl/mag   #mode field diameter @ 1/e**2

# fig = plt.figure(0)
# ax1=fig.add_subplot(3,1,1)
# ax2=fig.add_subplot(3,1,2)
# ax3=fig.add_subplot(3,1,3)

# ax1.imshow(slit[:,84:260])
# ax2.plot(wm, '-')
# ax3.plot(com-20, '-')

# ax1.set_ylabel('Position across \nthe slit ($\mu$m)')
# ax2.set_ylabel('Mode Field \nDiameter ($1/e^{2}, \mu$m)')
# #ax2.set_xlabel('Position along the slit ($\mu$m)')
# ax3.set_ylabel('Barycentre \n Position ($\mu$m)')
# ax3.set_xlabel('Position along the slit ($\mu$m)')

# plt.show()




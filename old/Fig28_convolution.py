from taskinit import *
from tasks import*
from matplotlib.pyplot import *
from cmath import sqrt
from operator import itemgetter

import os
import sys
import os.path
import math
import scipy.optimize
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *
import glob
import datetime
from astropy.units import beam


os.system('rm -rf *.ms *.im *.cl zero.* cont.* *.last *test* ft* Fig4* simobserve*')


##################################################### template image
ia.fromshape('zero.im',[512,512,1,1],overwrite=True)
cs=ia.coordsys()
cs.setunits(['rad','rad','','Hz'])
cell_rad=qa.convert(qa.quantity("0.1arcsec"),"rad")['value']
cs.setincrement([-cell_rad,cell_rad],'direction')
cs.setreferencevalue([qa.convert("03h32m39.53769s",'rad')['value'],qa.convert("-27d46m43.5673s",'rad')['value']],type="direction")
cs.setreferencevalue("228.1776275635GHz",'spectral')
cs.setincrement('15.626MHz','spectral') ## increased from 15.625 to avoid tclean errors for epoch=>10
ia.setcoordsys(cs.torecord())
ia.setbrightnessunit("Jy/pixel")
ia.done()


newim   =   ia.newimagefromimage(infile     =  'zero.im',
                                outfile    =  'cont.im')
data    =   newim.getchunk()
data[256,256,0,0]   =   1
newim.putchunk(data)
newim.close()


########################################## Imge convolution
exportfits(imagename='cont.im', fitsimage='cont.fits')
ia.open('cont.im')
ia.convolve2d(outfile = 'cont.conv.im', axes = [0,1], type = 'gauss', major = qa.quantity(0.488, 'arcsec'), minor = qa.quantity(0.463, 'arcsec'), pa=qa.quantity(-31.367, 'deg'))
ia.close()


##################################################### Simobserve
simobserve(skymodel     =   'cont.fits',
           inbright     =   '1Jy/pixel',
           indirection  =   'J2000 03h32m39.53769 -27d46m43.5673',
           incell       =   '0.1arcsec',
           incenter     =   '228.1776275635GHz',
           inwidth      =   '15.625MHz',
           project      =   'Fig4test',
           setpointings =   False,
           #
           ptgfile      =   'sim.ptg24.txt',
           #
           obsmode      =   'int',
           refdate      =   '2014/01/4',
           integration  =   '34.272s', #after scan avarage
           totaltime    =   '1', # number without units = time of repeating
           antennalist  =   "epoch1_hudf_spw2_scanav_field24.ms.cfg"
           )

os.system('cp -r Fig4test/Fig4test.epoch1_hudf_spw2_scanav_field24.ms.ms simobserve.ms')
tclean(vis = 'simobserve.ms', imagename = 'simobserve.ms', specmode =  'mfs', weighting     =   'natural',niter = 0, imsize = 512, cell = '0.1arcsec', pbcor = True)


###############################   ft and uvsub
exportfits(imagename='zero.im', fitsimage='zero.fits')
simobserve(skymodel     =   'cont.fits',
           inbright     =   '0Jy/pixel',
           indirection  =   'J2000 03h32m39.53769 -27d46m43.5673',
           incell       =   '0.1arcsec',
           incenter     =   '228.1776275635GHz',
           inwidth      =   '15.625MHz',
           project      =   'Fig4rep',
           setpointings =   False,
           #
           ptgfile      =   'sim.ptg24.txt',
           #
           obsmode      =   'int',
           refdate      =   '2014/01/4',
           integration  =   '34.272s', #after scan avarage
           totaltime    =   '1', # number without units = time of repeating
           antennalist  =   "epoch1_hudf_spw2_scanav_field24.ms.cfg"
           )

os.system('cp -r Fig4rep/Fig4rep.epoch1_hudf_spw2_scanav_field24.ms.ms ft_uvsub.ms')

os.system('rm -rf *.cl')
cl.addcomponent(dir='J2000 03h32m39.53769 -27d46m43.5673', flux=1, fluxunit='Jy',shape='point')
cl.rename('model.cl')
cl.done()

ft(vis='ft_uvsub.ms', spw='0:0~0',complist='model.cl',usescratch=True)
uvsub(vis='ft_uvsub.ms',reverse=True)
tclean(vis = 'ft_uvsub.ms', imagename = 'ft_uvsub.ms', specmode =  'mfs', weighting     =   'natural',niter = 0, imsize = 512, cell = '0.1arcsec', pbcor = True)


################################### flux distribusion
ia.open('cont.conv.im')
#ia.open('simobserve.ms.image')
#ia.open('ft_uvsub.ms.image')
data = ia.getchunk()
plt.ioff()
plt.clf()
plt.xlim([-0.2,1.2])
plt.xlabel('flux [Jy/B]', fontsize=18)
plt.ylabel('Npixels', fontsize=18)
plt.hist(data[~np.isnan(data)], bins = 30, alpha=0.2, normed=False, log=True, color='blue')
plt.legend(fontsize=18)
plt.savefig('test.pdf')
hist_flux = np.histogram(data[~np.isnan(data)], bins=30, range=None, normed=False, weights=None, density=None)
ia.done()

#b = data.ravel()
#test = np.sort(b)[::-1]

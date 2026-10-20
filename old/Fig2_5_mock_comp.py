# you need to install analysis utility as au
import glob
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
import matplotlib.gridspec as grids

Nbin = 50
def get_histogram(hist):
    array1 = np.insert(hist[0],0,0)
    array1 = np.append(array1,0)
    array2 = []
    bin_now = hist[1][1] - hist[1][0]
    for i, h in enumerate(hist[1][:-1]):
        array2.append(0.5*(hist[1][i+1] + hist[1][i]))
    array2 = np.insert(array2,0,array2[0]-bin_now)
    array2 = np.append(array2, array2[-1]+bin_now)
    stacked = np.vstack((array1,array2))
    return stacked

##############################  reset
os.system('rm -rf *dirtycube* sim* *rep* *last *listobs zero.im *.dat')
#del gain, replace, original

################################## original data
# Historical input: Field10_Spw0.ms from the archived ADF22 mock-data set.
vis_split   = 'Field10_Spw0.ms'
listobs(vis=vis_split, listfile=vis_split+'.listobs')

tb.open(vis_split, nomodify=False)
flagrow         =       tb.getcol('FLAG_ROW')
if False in flagrow:
    check       =       'available'
    flagging    =       True
else:
    check       =       'all flagged'
    flagging    =       False

name    = tb.colnames()
if 'CORRECTED_DATA' in name:
    column_now  =   'CORRECTED_DATA'
    datacolumn_now = 'corrected'
else:
    column_now  =   'DATA'
    datacolumn_now = 'data'

nchans      =   tb.getcol(column_now).shape[1]
original = tb.getcol(column_now)
gain = tb.getcol(column_now)
replace = tb.getcol(column_now)
tb.close()

'''
for chan in range(nchans):
    gain[0][chan][flagrow == 0]  = complex(np.nan, np.nan)
    gain[1][chan][flagrow == 0]  = complex(np.nan, np.nan)

gain.real[np.where(gain.real==0.)] = np.nan
gain.imag[np.where(gain.imag==0.)] = np.nan
'''

std_imag0 = np.nanstd(gain.imag[0], axis=(1))
mean_imag0 = np.nanmean(gain.imag[0], axis=(1))
std_real0 = np.nanstd(gain.real[0], axis=(1))
mean_real0 = np.nanmean(gain.real[0], axis=(1))
std_imag1 = np.nanstd(gain.imag[1], axis=(1))
mean_imag1 = np.nanmean(gain.imag[1], axis=(1))
std_real1 = np.nanstd(gain.real[1], axis=(1))
mean_real1 = np.nanmean(gain.real[1], axis=(1))


################################################# replaced data

size_now = np.size(flagrow)
for chan in range(nchans):
    #size_now = 2*np.size(gain.real[0][chan])
    replace.real[0][chan] = np.random.normal(loc = mean_real0[chan],scale =  np.sqrt(3)*std_real0[chan], size = size_now)
    replace.imag[0][chan] = np.random.normal(loc = mean_imag0[chan],scale =  np.sqrt(3)*std_imag0[chan], size = size_now)
    replace.real[1][chan] = np.random.normal(loc = mean_real1[chan],scale =  np.sqrt(3)*std_real1[chan], size = size_now)
    replace.imag[1][chan] = np.random.normal(loc = mean_imag1[chan],scale =  np.sqrt(3)*std_imag1[chan], size = size_now)
    '''
    j, k = 0, 0
    for i in range(size_now):
        if gain.real[0][chan][i] != 0.:
            replace.real[0][chan][i] = real_now[j]
            replace.real[1][chan][i] = real_now[j+1]
            j = j+2
        if gain.imag[0][chan][i] != 0.:
            replace.imag[0][chan][i] = imag_now[k]
            replace.imag[1][chan][i] = imag_now[k+1]
            k = k+2
    '''

vis_rep     = 'Field10_Spw0_rep.ms'
os.system('rm -rf vis_temp.ms')
os.system('cp -r '+vis_split+' vis_temp.ms')
tb.open('vis_temp.ms', nomodify=False)
tb.putcol(column_now, replace)
tb.flush
tb.close()
hanningsmooth(vis='vis_temp.ms', datacolumn=datacolumn_now, outputvis=vis_rep) # becareful about the datacolumn!!
listobs(vis=vis_rep, listfile=vis_rep+'.listobs')
os.system('rm -rf vis_temp.ms')


############################################### data obstaining and Imaging

uvnames = ['Field10_Spw0_rep.ms', 'Field10_Spw0.ms'] #glob.glob('*.ms')
for uvname in uvnames:
    tclean(vis=uvname,
            imagename=uvname+'.dirtycube',
            spw='0',
            specmode='cube',
            weighting='natural',
            datacolumn=datacolumn_now,
            width='15.625MHz',
            niter=0,
            cell = '0.1arcsec',
            imsize = 256,
            pblimit=0.8)

imnames = ['Field10_Spw0_rep.ms.dirtycube.image', 'Field10_Spw0.ms.dirtycube.image']# glob.glob('*.image')
imdata = []
chans, std, mean = [], [], []
for imname in imnames:
    ia.open(imname)
    datacube = ia.getchunk()
    datacube[datacube==0.0] = np.nan
    datacube *= 1000.
    imdata.append(datacube)
    chans.append(ia.shape()[3])
    std.append(np.nanstd(datacube, axis=(0,1,2)))
    mean.append(np.nanmean(datacube, axis=(0,1,2)))
    ia.close()

    result=imstat(imagename=imname)
    bmaj=imhead(imagename=imname, mode='list')['perplanebeams']['median area beam']['major']['value']
    bmin=imhead(imagename=imname, mode='list')['perplanebeams']['median area beam']['minor']['value']
    bpa=imhead(imagename=imname, mode='list')['perplanebeams']['median area beam']['positionangle']['value']
    print '## '+imname+': RMS = '+'{0:7.3f} mJy'.format(result['rms'][0]*1000.)
    print '##  Beamsize [arcsec]: {0:.3f} x {1:.3f} ({2:.3f}) ({2:.3f})'.format(bmaj,bmin,sqrt(bmaj*bmin),bpa)

uvdata = [replace.real, gain.real]
uvdata_i = [replace.imag, gain.imag]
histuv, histuv_i, histim = [], [], []
for gain, gain_i, cube in zip(uvdata, uvdata_i, imdata):
    hist = np.histogram(gain[~np.isnan(gain)], bins=Nbin, range=None, normed=False, weights=None, density=None)
    stacked = get_histogram(hist)
    histuv.append(stacked)

    hist = np.histogram(gain_i[~np.isnan(gain_i)], bins=Nbin, range=None, normed=False, weights=None, density=None)
    stacked = get_histogram(hist)
    histuv_i.append(stacked)


    hist = np.histogram(cube[~np.isnan(cube)], bins=Nbin, range=None, normed=False, weights=None, density=None)
    stacked = get_histogram(hist)
    histim.append(stacked)

hist = np.histogram(original[~np.isnan(original)], bins=Nbin, range=None, normed=False, weights=None, density=None)
stacked = get_histogram(hist)


################################################# Figure

plt.ioff()
plt.clf()
fig = plt.figure(figsize=(10, 3)) # inch
gs = grids.GridSpec(12,11) # 0,0 => left top

pl = plt.subplot(gs[:,0:3])
#pl.text(2, 0.9, 'Visibility Plane', size = 10, color ='black')
pl.set_xlim([-26.,26])
pl.set_ylim([0.8,1e8])
pl.set_xlabel('flux [Jy]')
pl.set_ylabel('Npixels')
#pl.set_xscale('log')
pl.set_yscale('log')
pl.step(stacked[1], stacked[0], where='mid', linewidth=1, color='black',zorder=-1, label='Replaced')
pl.step(histuv[0][1], histuv[0][0], where='mid', linewidth=1, color='blue',zorder=-1, label='Replaced')
pl.step(histuv[1][1], histuv[1][0], where='mid', linewidth=1, color='red',zorder=-1, label='Original')
#pl.step(histuv_i[0][1], histuv_i[0][0], where='mid', linewidth=1, color='blue', linestyle='dashed', zorder=-1, label='Replaced')
#pl.step(histuv_i[1][1], histuv_i[1][0], where='mid', linewidth=1, color='red', linestyle='dashed', zorder=-1, label='Original')
#pl.legend(bbox_to_anchor=(1,0.8), loc='upper right')


pm = plt.subplot(gs[:,4:7])
#pm.text(2, 0.9, 'Image Plane', size = 10, color ='black')
pm.set_xlim([-8.,8.])
pm.set_ylim([0.8,1e6])
pm.set_xlabel('flux [mJy/B]')
pm.set_ylabel('Npixels')
#pm.set_xscale('log')
pm.set_yscale('log')
pm.step(histim[0][1], histim[0][0], where='mid', linewidth=1, color='blue',zorder=-1, label='Replaced')
pm.step(histim[1][1], histim[1][0], where='mid', linewidth=1, color='red',zorder=-1, label='Original')
#pm.legend(bbox_to_anchor=(1,0.8), loc='upper right')

pr = plt.subplot(gs[:,8:11])
#pr.text(2, 0.9, 'Image Plane', size = 19, color ='black')
pr.set_ylim([0.6,1.])
pr.set_xlim([2,120])
pr.set_xlabel('frequency channels')
pr.set_ylabel('RMS value [mJy]')
pr.plot(np.arange(127), std[0], color='blue', label='Replaced')
pr.plot(np.arange(127), std[1], color='red', label='Original')
#pr.legend(bbox_to_anchor=(1,0.8), loc='upper right')

plt.savefig('Fig2_5_compMockdata.pdf', bbox_inches='tight', pad_inches=0)

'''

 for chan_now in xrange(0,chans):
    if std_data[chan_now] > 0.0:
        data[:,:,0,chan_now] = data[:,:,0,chan_now]/std_data[chan_now]
    else:
        data[:,:,0,chan_now] = np.nan


plt.ylim([1,1e8])
plt.xlim([-8,8])
plt.xlabel('S/N ratio [sigma]')
plt.ylabel('Npixels')
plt.hist(data[~np.isnan(data)], bins = Nbin, alpha=0.2, normed=False, log=True, color='green',  histtype="stepfilled", label=name)
'''

'''
#################################################  Simobserve
ms.open(vis_split)
scaninfo = ms.getscansummary()
scanid = scaninfo.keys()
inttime = scaninfo[str(scanid[0])]['0']['IntegrationTime']
time_begin = qa.time(qa.quantity(scaninfo[str(scanid[0])]['0']['BeginTime'],'d'), form='ymd')
date_now, hour = time_begin[0].rsplit("/",1)[0], time_begin[0].rsplit("/",1)[1]
ms.close()

field_now = 10
Spw = 0
print str(field_now)+' '+str(Spw)+' '+str(mean_real)+' '+str(std_real)+' '+str(mean_imag)+' '+str(std_imag)+' '+str(inttime)+' '+str(date_now)+' '+str(nchans)+' '+str(check)+' '+str(flagging)+' '+column_now

# Historical input: SB0_Field10_Spw0.txt from the archived pointing files.
filename    =   'pointingfiles/SB0_Field10_Spw0.txt'
textfile    =   open(filename, "r")
line_start  =   1
pointings   =   textfile.readlines()[line_start:]
RA, Dec     =   [], []
id_field    =   [10]
for pointing in pointings:
    temp = pointing.replace('\n','').split()
    RA.append(temp[1])
    Dec.append(temp[2])
textfile.close()

# Historical input: the corresponding archived array-configuration files.
nu_center   =   '254.0000000GHz'
RA          =   RA[0].split(':')
Dec         =   Dec[0].split('.')
num_now     =   10
date_now    =   '2014/06/06'
hourang     =   float(num_now*30/3600.)
del_nu      =   '15.625MHz'#128chans
nchans      =   128
inttime     =   40

#ms.open(vis_split)
#scaninfo = ms.getscansummary()
#scanid = scaninfo.keys()
#inttime = scaninfo[str(scanid[0])]['0']['IntegrationTime']
#ms.close()

ia.fromshape('zero.im',[512,512,1,nchans],overwrite=True)
cs=ia.coordsys()
cs.setunits(['rad','rad','','Hz'])
cell_rad=qa.convert(qa.quantity('0.1arcsec'),"rad")['value']
cs.setincrement([-cell_rad,cell_rad],'direction')
cs.setreferencevalue([qa.convert(RA[0]+'h'+RA[1]+'m'+RA[2]+'s','rad')['value'],qa.convert(Dec[0]+'d'+Dec[1]+'m'+Dec[2]+'s','rad')['value']],type="direction")
cs.setreferencevalue(nu_center,'spectral')
cs.setincrement(del_nu,'spectral')
ia.setcoordsys(cs.torecord())
ia.setbrightnessunit("Jy/pixel")
ia.done()

simobserve( skymodel    =   'zero.im',
            project     =   'sim',
            setpointings=   False,
            ptgfile     =   'pointingfiles/SB0_Field10_Spw0.txt',
            refdate     =   date_now,
            obsmode     =   'int',
            integration =   str(inttime)+'s',
            totaltime   =   '9',
            hourangle   =   str(hourang)+'h',
            antennalist =   'configurationfiles/SB0_Field10_Spw0.cfg',
            overwrite   =   True)

## add noise to noise-free ms
os.system('rm -rf sim/*.noisy.ms')
rms_vis = np.abs(sqrt(float(std_real)*float(std_imag)))
sm.openfromms('sim/sim.SB0_Field10_Spw0.ms')
sm.setseed(num_now)
sm.setnoise(mode='simplenoise',simplenoise='%fJy' % rms_vis)
sm.corrupt()
sm.done()
os.system('mv sim/sim.SB0_Field10_Spw0.ms sim.SB0_Field10_Spw0.ms')
'''


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#                                                                                                                #
#  WHAT NEED TO BE INSTALLED - CASA : https://casa.nrao.edu/casa_obtaining.shtml                                #
#                                                                                                                #
#  HOW TO RUN THIS SCRIPT ---- 1. set your parameters in this file    (line 48~)                                    #
#                              2. In CASA commandline, 'cd' to the same dir. to this .py file                    #
#                                      (example) CASA<2> : cd path_to_scriptfile_directry                        #
#                              3. Run this script                                                                #
#                                      (example) CASA<3> : execfile('linesearch_170222.py')                        #
#                              4. Input field name                                                                #
#                              5. Check output file linesearch.log                                                    #
#                                                                                                                #
#  TEST RUN ------------------- Please find test run data files from:                                            #
#                                      https://... ... ...                                                        #
#                                                                                                                #
#  MAIN REFERENCES :                                                                                            #
#    N. H. Hayatsu and ADF22 team. 2017, accepted to PASJ,   (method for line search)                            #
#    CASA Toolkit Reference Manual (web page):                                                                    #
#              http://www.eso.org/~scastro/ALMA/casa/active/CasaRef/CasaRef.html#GeneralPackage.html            #
#                                                                                                                #
#  HISTORY :                                                                                                    #
#   **/07/2017  version 0.2   (varid for CASA ver.4.7.0) Natsuki Hayatsu                                        #
#            - include completeness chack routine                                                                #
#   16/06/2017  version 0.12  (varid for CASA ver.4.7.0) Natsuki Hayatsu                                        #
#       updates from previous version                                                                            #
#           - reduce input parameters; calc. using msmd and ma tools                                            #
#                estimated survey area ref: https://almascience.eso.org/about-alma/alma-basics                    #
#           - clean -> make use of task tclean                                                                    #
#           - uvcontsub using excludechans option                                                                #
#   26/05/2017  version 0.11 (runs on CASA ver.4.*.* and 5.0.*) Kanako Sugimoto (modification log for Hayatsu-san)    #
#       updates from previous version                                                                            #
#           - introduced smothing using scipy.signal.convolve to speed up processing                            #
#           - fixed hadling of right edge channels in smoothing                                                    #
#           - explicitly deleted cubes on memory that are not used in following stages to save memory consumption    #
#           - robust Measurementset summary using msmd tool                                                        #
#   17/03/2017  version 0.0 ( varid for CASA ver.4.*.*)                                                            #
#       updates from N.H.Hayatsu et al. (2017) accepted to PASJ                                                    #
#           - search region : rectangle -> high sensitivity area (users can define the threshold)                #
#           - SNcube        : normarized by RMS -> by standard deviation                                        #
#           - Clumpfind     : original algorism                                                                    #
#                                                                                                                #
#                                         Natsuki H. Hayatsu                                                     #
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

from taskinit import *
from tasks import*

import os
import sys
import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import scipy.signal


## input parameters ##
dir_name    =    '../SXDF_band8/concat_all.ms'
low            =    4.0                #lower limit of S/N ratio for search
vel_max        =    400                # maximum smoothing palameter in km/s
vel_min        =    100                # 0 = original resolution
Nbin_vel    =    50                # number of smoothing data

MASK        =    0.8
CHANS        =    0.9
BUF_XPIXEL    =    10
BUF_YPIXEL    =    10
V_RESO        =    1                # minimum velocity resolution  km/s
NPIXMIN        =    5
NCLUMP_TOT    =    9999
NSOURCE_TOT    =    9999
NCORE        =    1
ANTENNA        =    12e-3            # 12m -> 12e-3km/s

## parameters controlling output ##
OUTPUT_SMCUBE        =    1       #smoothing image cube
OUTPUT_SNCUBE        =    1
OUTPUT_ARTIFICIAL    =    1        #to be updated

#if you want skip some steps: 0(don't skip) or 1(skip)
STEP123        =    np.array([0, 0, 1], dtype = 'i1')
STEP4567    =    1
STEP8        =    1

#-------------------------------------------------------------------Step 0  Read obsearvational/linesarch parameters
thesteps = []
step_title = {0: ' Read obsearvational/linesarch parameters',
              1: ' Generate contnuum image data',
              2: ' Generate contnuum subtracted image data cube',
              3: ' Masking',
              4: ' Spectral smoothing',
              5: ' ClumpFind in max S/Ncube',
              6: ' Comtamination check',
              7: ' Line Search Summary',
              8: ' Completeness check'}

mystep = 0
casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
print '#Step', mystep, step_title[mystep]

ms_dir0 = str(dir_name)
ms_dir1  = str(dir_name)+'.contsub'

print '#                   ms_dir0 : %s' % ms_dir0

sample = 'linesearch/param/listobs.cat'

input_ms = ms_dir0

if not os.path.exists(input_ms):
    if os.path.exists(ms_dir1):
        input_ms = ms_dir1
    else:
        raise Exception, "Could not find an input MS"

if not os.path.exists(sample):
    os.system('mkdir linesearch;mkdir linesearch/param')
    listobs(vis=input_ms,listfile=sample)

msmd.open(input_ms)
im.open(input_ms)
try:

    # input field name
    field_names = msmd.fieldsforintent('*', asnames=True)
    print '#                   field names:'
    for name in field_names:
        print '#                        %s' % name
    print '\n#                    Input char. field_name ...'
    field_name=raw_input()


    print '#Step', mystep, step_title[mystep]
    print '#  user parameters   ms_dir0 : %s' % ms_dir0
    print '#                    ms_dir1 : %s' % ms_dir1
    print '#                    field name : %s' % field_name
    print '#                    search  > %.1f sigma clumps' % low
    print '#                    line search  %d ~ %d km/s, with %d bins' % (vel_min, vel_max, Nbin_vel)
    print '#                        (0 = original resolution)'


    # spw
    targetspw = msmd.spwsforintent('*TARGET*')
    num_spw = len(targetspw)
    print '\n#                    number of spectral window : %d' % num_spw

    # observed frequency
    nu_list = [ msmd.chanfreqs(spwid, 'GHz')[0] for spwid in xrange(num_spw) ]
    nu_center = (max(nu_list)+min(nu_list))/2.0
    print '#                    nu_center : %f GHz (%f ~ %f GHz)' % (nu_center, min(nu_list), max(nu_list))

    # original velocity resolution
    vel_ori = [ msmd.chanres(spwid, 'km/s', asvel=True) for spwid in xrange(num_spw) ]
    for spw_now in range(num_spw):
        print '#                    original velocity resolusion : %f km/s' % (np.mean(vel_ori[spw_now]))

    # output channels
    chans_array = [msmd.nchan(spw_now) for spw_now in xrange(num_spw)]
    for spw_now in range (num_spw):
        print '#                    chans spw%d: %d' % (spw_now, chans_array[spw_now])

    # obtain intent field IDs
    field_IDs = msmd.fieldsforname(field_name)
    num_fields = len(field_IDs)
    print '#                    number of fields : %d' % num_fields
    print '#                    field ID:  '
    for ID in field_IDs:
        print '#                        %d' % ID


	# cell size
    	cell = im.advise()[2]['value']    # beam_size / 2 [arcsec] <--------------------------check it ! ~ / 3 ?
    	print '#                    cell size: %.2f arcsec ~ beam size/2' % cell

	# estimated survey area
	fov_singlepoint =	76 / ANTENNA / nu_center / 60	# arcmin
	print '#                    fov for single beam: %.2f arcsec' % (fov_singlepoint*60)

	dir_antenna = [msmd.phasecenter(ID) for ID in field_IDs]
	dir_val_x = [msmd.phasecenter(ID)['m1']['value'] for ID in field_IDs]	# radian
	dir_val_y = [msmd.phasecenter(ID)['m0']['value'] for ID in field_IDs]	# radian

	id_a = np.where(np.array(dir_val_x) == min(dir_val_x))
	id_b = np.where(np.array(dir_val_y) == max(dir_val_y))
	id_c = np.where(np.array(dir_val_x) == max(dir_val_x))

	dir_a = dir_antenna[id_a[0][0]]
	dir_b = dir_antenna[id_b[0][0]]
	dir_c = dir_antenna[id_c[0][0]]

	h_sep = me.separation(dir_b, dir_c)['value']*60	# arcmin
	v_sep = me.separation(dir_a, dir_b)['value']*60 # arcmin

	if h_sep < 1e-10: h_sep = 0.0
	if v_sep < 1e-10: v_sep = 0.0

	phi = math.radians(me.posangle(dir_b, dir_c)['value'])
	theta = math.radians(me.posangle(dir_b, dir_a)['value'])
	print '#                    position angle: %.2f rad and %.2f rad' % (phi, theta)

	if h_sep == 0 and v_sep == 0:
		survey_area = fov_singlepoint * fov_singlepoint
		gridder		=	'standard'
		#deconvolver=	'multiscale'
		deconvolver=	'hogbom'#!-----------------------------------170815
		h_sep = fov_singlepoint
		v_sep = fov_singlepoint
		imsize_x = (int) (fov_singlepoint*60 / cell + 100)
		imsize_y = (int) (fov_singlepoint*60 / cell) + 100
		print '#                    survey area(single point): %.2f ^2 = %.2f arcmin^2' \
									% (fov_singlepoint, survey_area)

	else:
		survey_area = (v_sep + fov_singlepoint) * (h_sep + fov_singlepoint)
		gridder		=	'mosaic'
		deconvolver=	'hogbom'
		imsize_x = int((h_sep*math.fabs(math.sin(phi)) + v_sep*math.fabs(math.sin(theta)))*60/cell + 100)
		imsize_y = int((h_sep*math.fabs(math.cos(phi)) + v_sep*math.fabs(math.cos(theta)))*60/cell + 100)
		print '#                    survey area: %.2f x %.2f = %.2f arcmin^2' \
					% (h_sep+fov_singlepoint, v_sep+fov_singlepoint, survey_area)


	print '#                    imsize: %d x %d' % (imsize_x, imsize_y)



    # estimated flux sensitivity
    im.selectvis(field = field_name, spw= '')
    im.weight(type='natural')
    #sensitivity=im.apparentsens()
    #sigma_cont = sensitivity[1]*1e3* math.sqrt(num_field) # mJy/beam
    print '#                    sigma_cont: %.3f mJy/beam' % sigma_cont

finally:
    msmd.close()
    im.close()

print ''

#-------------------------------------------------------------------Step 1  Generate contnuum image data

if not os.path.exists('linesearch/original'):
    os.system('mkdir linesearch/original')

mystep = 1
if STEP123[0] == 0:
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print ('#Step', mystep, step_title[mystep])


    os.system('rm -rf linesearch/original/clean.cont*')

    tclean(vis= ms_dir0,
    imagename = 'linesearch/original/clean.cont',
    field = field_name,
    gridder=gridder,
    specmode='mfs',
    spw = '',
    cell = '%.2f arcsec' % cell,
    imsize=[imsize_x, imsize_y],
    interpolation = 'nearest',
    weighting='natural',
    niter=0,
    ##threshold='%.3f mJy' % sigma_cont,
    pbcor=True,
    parallel = False,    #to be updated
    interactive=False)

else:
    print '#Step', mystep, step_title[mystep], ' is skipped.'


print ''


#-------------------------------------------------------------------Step 2  Generate contnuum subtracted datacube
mystep = 2
if(STEP123[1] == 0):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#Step', mystep, step_title[mystep]


    if not os.path.exists('linesearch/original'):
        os.system('mkdir linesearch/original')

    os.system('rm -rf linesearch/original/spw*.cube.contsub*')

    for spw_now in range(num_spw):

        if V_RESO > np.mean(vel_ori[spw_now]):
            V_RESO = np.mean(vel_ori[spw_now])

        os.system('rm -rf '+ms_dir1)

        uvcontsub(vis = ms_dir0,
            field = field_name,
            fitspw = '%d' % spw_now,
            excludechans = True,
            combine = 'spw',
            solint = 'int',
            fitorder = 0,
            spw = '%d' % spw_now,
            want_cont = False)

        tclean(vis=ms_dir1,
            imagename = 'linesearch/original/spw%d.cube.contsub' % spw_now,
            field = field_name,
            specmode='cube',
            spw = '0',
            interpolation = 'nearest',
            gridder=gridder,
            chanchunks = 4, #! -----------------------
            imsize=[imsize_x, imsize_y],
            cell = '%.2f arcsec' % cell,
            #width='%f km/s' % V_RESO,
            niter=0,
            weighting='natural',
            pbcor=True,
            outframe='BARY',
            #threshold='0.01mJy',
            parallel = False,
            interactive = False)

else:
    print '#Step', mystep, step_title[mystep], ' is skipped.'


print ''


#-------------------------------------------------------------------Step 3  Masking

mystep = 3
if(STEP123[2] == 0):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#Step', mystep, step_title[mystep]


    if os.path.exists('linesearch/masked'):
        os.system('rm -rf linesearch/masked')
    os.system('mkdir linesearch/masked')

    for spw_now in range(num_spw):
        immath(imagename='linesearch/original/spw%d.cube.contsub.image' % spw_now,
               expr = 'IM0*1000',
               #mask= '"linesearch/original/spw%d.cube.contsub.pb" >= %.1f' % (spw_now,MASK),
               outfile='linesearch/masked/spw%d.cube.contsub.image' % spw_now)

        #makemask(mode='copy',
        #         'linesearch/masked/spw%d.cube.contsub.image' % spw_now,
        #         inpmask='"linesearch/original/spw%d.cube.contsub.pb" >= %.1f' % (spw_now,MASK),
        #         outfile='linesearch/masked/spw%d.cube.contsub_%.1f.image' % (spw_now, MASK))


else:
    print '#Step', mystep, step_title[mystep], ' is skipped.'


print ''


#------------------------------------------------------------------- Step 4  Spectral smoothing
for spw_now in range(num_spw): #roop1 ... for spectral window
    if STEP4567 == 1:
        print '#Step 4 -- 5 is skipped.'
        print ''

        break

    mystep = 4
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#roop1 for spw. : spw%d.' % spw_now
    print ''
    print '#    Step', mystep, step_title[mystep]


    if os.path.exists('linesearch/masked/spw%d.cube.contsub.masked.image' % spw_now):
        masked_file = 'linesearch/masked/spw%d.cube.contsub.masked.image' % spw_now
    else:
        masked_file = 'linesearch/masked/spw%d.cube.contsub_%.1f.image' % (spw_now, MASK)

	stat = imstat(imagename=masked_file,axes=[0,1],
	logfile='linesearch/masked/spw%d.cube.contsub.masked.imstat' % spw_now)


    stat_sigma = stat["sigma"]
    max_sigma = np.nanmax(stat_sigma)

    chans = chans_array[spw_now]
    start_chan = chans
    sm_chans = 0
    end_chan = 0
    bad_slice = np.zeros([chans], dtype = 'int')

    for chan_now in xrange(0,chans):
        if  stat["sigma"][chan_now] <= max_sigma * CHANS and stat["sigma"][chan_now] > 0.0:
            if start_chan > chan_now : start_chan = chan_now
            if end_chan < chan_now :   end_chan   = chan_now
        else:
            bad_slice[chan_now] = -1 # not used for now

    sm_chans = end_chan - start_chan + 1

    print '#        roop2 for smoothing.'

    #if V_RESO > np.mean(vel_ori[spw_now]):
    V_RESO = np.mean(vel_ori[spw_now])

    sm_min = int(vel_min/V_RESO)
    sm_max = int(vel_max/V_RESO)
    d_sm   = int(Nbin_vel/V_RESO)


    for sm_now in range(sm_min,sm_max+d_sm,d_sm): #roop2 ... for smoothing parameter
        chan_first = int(start_chan +  sm_now/2)
        chan_last  = int(end_chan - sm_now/2)
        print '#        spw%d sm_now = %d chans=[%d, %d]' % (spw_now, sm_now, chan_first, chan_last)


        smoothed = np.zeros([imsize_x,imsize_y,chans], dtype = 'float')
        SNcube = np.zeros([imsize_x,imsize_y,chans], dtype = 'float')
        std = np.zeros([chans], dtype = 'float')

        original = imval(masked_file,
                 chans='%d~%d' % (start_chan, end_chan),
                 box='0,0,%d,%d' % (imsize_x-1, imsize_y-1))
        cube = original['data']
        del original

        #KS NOTE: the following lines assume frequency axis to be the 3rd one.
        nvalid_chan = end_chan-sm_now-start_chan+1 # the number of channels used in smoothed cube
        if sm_now == 0: # NO SMOOTHING
            smoothed[:,:,start_chan:end_chan+1] = cube[:,:,0:sm_chans]
        else:
            kernel = np.ones((1,1,sm_now), dtype=np.float64)/float(sm_now)
            smoothed[:,:,chan_first:nvalid_chan+chan_first] = scipy.signal.convolve(cube, kernel, mode='same')[:,:,sm_now/2:nvalid_chan+sm_now/2]
        del cube
        std = np.std(smoothed, axis=(0,1))
        for chan_center in xrange(chan_first, chan_first+nvalid_chan):
                        SNcube[:,:,chan_center] = smoothed[:,:,chan_center] / std[chan_center]


        if OUTPUT_SMCUBE == 0:
            if spw_now == 0 and sm_now == sm_min:
                if os.path.exists('linesearch/smoothed'):
                    os.system('rm -rf linesearch/smoothed')
                os.system('mkdir linesearch/smoothed')

            #try
            ia.fromarray(pixels=smoothed,
                outfile='linesearch/smoothed/spw%d.cube.contsub.masked.sm%d.image' % (spw_now, sm_now))
            #finally
            # if ia.isopen ...
            ia.done()

        print ''

        del smoothed

#-------------------------------------------------------------------Step 5  Clumpfind in S/Ncube
        if sm_now == sm_min:
            max_SNcube = SNcube
            inv_SNcube = SNcube
        else:
            max_SNcube = np.fmax(SNcube, max_SNcube)#matching in 3D
            inv_SNcube = np.fmin(SNcube, inv_SNcube)#inverted max SNcube before x -1

            #sum up and mean

        if spw_now == 0 and sm_now == sm_min:
            if os.path.exists('linesearch/SNcube'):
                os.system('rm -rf linesearch/SNcube')
            os.system('mkdir linesearch/SNcube')

        if OUTPUT_SNCUBE == 0:
            print '#    Generate S/N cube ....'
            print ''

            ia.fromarray(pixels=SNcube,
                outfile='linesearch/SNcube/spw%d.SNcube.contsub.masked.sm%d.image' % (spw_now, sm_now))
            ia.done()


        if spw_now == 0 and sm_now == sm_min:
            if os.path.exists('linesearch/Fig'):
                os.system('rm -rf linesearch/Fig')
            os.system('mkdir linesearch/Fig')

        # plot SNdist
        SNcube_1dim =  SNcube[:,:,chan_first:chan_last].reshape((imsize_x*imsize_y*(chan_last-chan_first)))
        #if spw_now == 0 and sm_now == sm_min:
        #    plot_ymax = np.size(np.where(SNcube_1dim == 0))
        #    plot_ymax = 1000*(int(plot_ymax/1000) + 5)

        plt.ioff()
        plt.clf()
        #plt.ylim([0.5,plot_ymax])
        plt.ylim([0.5,1e6])
        plt.xlim([-8,8])
        plt.xlabel('S/N ratio [sigma]')
        plt.ylabel('Npixels')
        plt.hist(SNcube_1dim[~np.isnan(SNcube_1dim)], bins = 100, log = True)
        figname = 'linesearch/Fig/SNdist.spw%d.sm%d.png' % (spw_now, sm_now)
        plt.savefig(figname)
        del SNcube, SNcube_1dim

        # end  for roop 2

    mystep = 5
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#    Step', mystep, step_title[mystep]


    print ''

    max_SNcube[max_SNcube < 0.0] = 0.0
    inv_SNcube[inv_SNcube > 0.0] = 0.0
    max_SNcube += inv_SNcube
    del inv_SNcube

    SNcube_1dim = max_SNcube[:,:,chan_first:chan_last].reshape((imsize_x*imsize_y*(chan_last-chan_first)))
    plt.clf()
    #plt.ylim([0.5,plot_ymax])
    plt.ylim([0.5,1e6])
    plt.xlim([-8,8])
    plt.xlabel('S/N ratio [sigma]')
    plt.ylabel('Npixels')
    figname = 'linesearch/Fig/SNdist.spw%d.max.png' % spw_now
    plt.hist(SNcube_1dim[~np.isnan(SNcube_1dim)], bins = 100, log = True, )
    plt.savefig(figname)
    del SNcube_1dim

    ia.fromarray(pixels=max_SNcube,
        outfile='linesearch/SNcube/spw%d.maxSNcube.contsub.masked.image' % spw_now)
    ia.done()


    for chan_now in range(start_chan+1,end_chan+1):#matching in 2D
        if chan_now == start_chan + 1 :
            max_SNmap_xy = np.fmax(max_SNcube[:,:,chan_now-1], max_SNcube[:,:,chan_now])
            min_SNmap_xy = np.fmin(max_SNcube[:,:,chan_now-1], max_SNcube[:,:,chan_now])

        else:
            max_SNmap_xy = np.fmax(max_SNmap_xy, max_SNcube[:,:,chan_now])
            min_SNmap_xy = np.fmin(min_SNmap_xy, max_SNcube[:,:,chan_now])

    for i in range(1,imsize_x):
        if i == 1 :
            max_SNmap_yz = np.fmax(max_SNcube[i-1,:,:], max_SNcube[i,:,:])
        else:
            max_SNmap_yz = np.fmax(max_SNmap_yz, max_SNcube[i,:,:])

    for j in range(1,imsize_y):
        if j == 1 :
            max_SNmap_xz = np.fmax(max_SNcube[:,j-1,:], max_SNcube[:,j,:])
        else:
            max_SNmap_xz = np.fmax(max_SNmap_xz, max_SNcube[:,j,:])


    plt.clf()
    plt.imshow(max_SNmap_xy, cmap=plt.cm.rainbow, vmin=-6.0, vmax=6.0)#!  vmin, vmax
    plt.xlabel('y pixel')
    plt.ylabel('x pixel')
    plt.colorbar().set_label('S/N ratio')
    figname = 'linesearch/Fig/max_SNmap_xy_%d.png' % spw_now
    plt.savefig(figname)

    plt.clf()
    plt.imshow(min_SNmap_xy, cmap=plt.cm.rainbow, vmin=-6.0, vmax=6.0)#!  vmin, vmax
    plt.xlabel('y pixel')
    plt.ylabel('x pixel')
    plt.colorbar().set_label('S/N ratio')
    figname = 'linesearch/Fig/min_SNmap_xy_%d.png' % spw_now
    plt.savefig(figname)

    plt.clf()
    plt.imshow(max_SNmap_xz, cmap=plt.cm.rainbow, vmin=-6.0, vmax=6.0)
    plt.xlabel('channels')
    plt.ylabel('x pixel')
    plt.colorbar().set_label('S/N ratio')
    figname = 'linesearch/Fig/max_SNmap_xz_%d.png' % spw_now
    plt.savefig(figname)

    plt.clf()
    plt.imshow(max_SNmap_yz, cmap=plt.cm.rainbow, vmin=-6.0, vmax=6.0)
    plt.xlabel('channels')
    plt.ylabel('y pixel')
    plt.colorbar().set_label('S/N ratio')
    figname = 'linesearch/Fig/max_SNmap_yz_%d.png' % spw_now
    plt.savefig(figname)


    print '#    search clump with  >= %.1f sigma '% low


    max_SNcube[(-1.0*low  < max_SNcube) & (max_SNcube < low)] = 0.0
    max_SNmap_xy[(-1.0*low  < max_SNmap_xy) & (max_SNmap_xy < low)] = 0.0

    id_SNcube = np.zeros([imsize_x,imsize_y,chans], dtype = 'int')
    id_SNmap = np.zeros([imsize_x,imsize_y], dtype = 'int')

    Nclump = Nsource = Nclump_inv = Nsource_inv = 0

    #positive
    for j in xrange(1,imsize_y):
        for i in xrange(0,imsize_x-1):
            if max_SNmap_xy[i,j] >= low: # when meet new clump at (i+k, j)
                if id_SNmap[i-1,j-1] > 0: id_SNmap[i,j] = id_SNmap[i-1,j-1]
                elif id_SNmap[i-1,j] > 0: id_SNmap[i,j] = id_SNmap[i-1,j]
                elif id_SNmap[i,j-1] > 0: id_SNmap[i,j] = id_SNmap[i,j-1]
                elif id_SNmap[i+1,j-1] > 0: id_SNmap[i,j] = id_SNmap[i+1,j-1]
                else:
                    Nsource += 1
                    id_SNmap[i, j] = Nsource

    temp_id = 0
    for source in range(1,Nsource+1):
        if np.size(np.where(id_SNmap == source))/2 >= NPIXMIN:
            temp_id += 1
            id_SNmap[id_SNmap == source] = temp_id
        else:
            id_SNmap[id_SNmap == source] = 0

    Nsource = temp_id

    temp_id = 0
    for source in range(1, Nsource+1):
        where_source = np.array(np.where(id_SNmap == source))
        y_len = len(where_source[0])
        for k in range(start_chan+1, end_chan+1):
            for y_now in range(0,y_len):
                i_now = where_source[0,y_now]
                j_now = where_source[1,y_now]
                if max_SNcube[i_now,j_now,k] >= low :
                    if id_SNcube[i_now,j_now,k-1] > 0 :
                        id_SNcube[i_now,j_now,k] = id_SNcube[i_now,j_now,k-1]
                    else:
                        temp_id += 1
                        id_SNcube[i_now, j_now,k] = temp_id

    clump_posi = 0
    for temp_clump in range(1,temp_id+1):
        if np.size(np.where(id_SNcube == temp_clump))/2 >= NPIXMIN:
            clump_posi += 1
            id_SNcube[id_SNcube == temp_clump] = clump_posi
        else: id_SNcube[np.where(id_SNcube == temp_clump)] = 0

    Nclump = clump_posi

    #negative
    for j in xrange(1,imsize_y):
        for i in xrange(0,imsize_x-1):
            if min_SNmap_xy[i,j] <= -1.0*low: # when meet new clump at (i+k, j)
                if id_SNmap[i-1,j-1] < 0: id_SNmap[i,j] = id_SNmap[i-1,j-1]
                elif id_SNmap[i-1,j] < 0: id_SNmap[i,j] = id_SNmap[i-1,j]
                elif id_SNmap[i,j-1] < 0: id_SNmap[i,j] = id_SNmap[i,j-1]
                elif id_SNmap[i+1,j-1] < 0: id_SNmap[i,j] = id_SNmap[i+1,j-1]
                else:
                    Nsource_inv += 1
                    id_SNmap[i, j] = -1*Nsource_inv

    temp_id = 0
    for source in range(1,Nsource_inv+1):
        if np.size(np.where(id_SNmap == -1*source))/2 >= NPIXMIN:
            temp_id += 1
            id_SNmap[id_SNmap == -1*source] = -1*temp_id
        else:
            id_SNmap[id_SNmap == -1*source] = 0

    Nsource_inv = temp_id

    temp_id = 0
    for source in range(1, Nsource_inv+1):
        where_source = np.array(np.where(id_SNmap == -1*source))
        y_len = len(where_source[0])
        for k in range(start_chan+1, end_chan+1):
            for y_now in range(0,y_len):
                i_now = where_source[0,y_now]
                j_now = where_source[1,y_now]
                if max_SNcube[i_now,j_now,k] <= -1.0*low :
                    if id_SNcube[i_now,j_now,k-1] < 0 :
                        id_SNcube[i_now,j_now,k] = id_SNcube[i_now,j_now,k-1]
                    else:
                        temp_id += 1
                        id_SNcube[i_now, j_now,k] = -1*temp_id

    clump_nega = 0
    for temp_clump in range(1,temp_id+1):
        if np.size(np.where(id_SNcube == -1*temp_clump))/2 >= NPIXMIN:
            clump_nega += 1
            id_SNcube[id_SNcube == -1*temp_clump] = -1*clump_nega
        else: id_SNcube[np.where(id_SNcube == -1*temp_clump)] = 0

    Nclump_inv = clump_nega

    print '#    Nsource %d Nsource_inv %d' % (Nsource, Nsource_inv)
    print '#    Nclump %d Nclump_inv %d' % (Nclump, Nclump_inv)



    # plot SN map id #
    plt.clf()
    plt.imshow(id_SNmap, cmap=plt.cm.rainbow, vmin=-1*Nsource_inv, vmax=Nsource)
    plt.colorbar().set_label('id')
    plt.xlabel('y pixel')
    plt.ylabel('x pixel')
    figname ='linesearch/Fig/id_SNmap_%d.png' % spw_now
    plt.savefig(figname)

    # save SN cube id #
    if spw_now == 0:
        if os.path.exists('linesearch/clumpfind'):
            os.system('rm -rf linesearch/clumpfind')
        os.system('mkdir linesearch/clumpfind')

    ia.fromarray(pixels=id_SNcube,
        outfile='linesearch/clumpfind/spw%d.cube.id.image' % spw_now)
    ia.done()


    clfind = np.zeros([Nclump,8], dtype = 'float')
    clfind_inv = np.zeros([Nclump_inv,6], dtype = 'float')

    if spw_now == 0:
        clfind_res = np.zeros([NCLUMP_TOT,8,num_spw], dtype = 'float')
        clfind_res_inv = np.zeros([NCLUMP_TOT,6,num_spw], dtype = 'float')

    print '#    read pbcor data...'


    pbcor = imval('linesearch/original/spw%d.cube.contsub.image.pbcor' % spw_now, box='0,0,%d,%d' % (imsize_x-1,imsize_y-1))
    data_flux = pbcor['data']*1e3 # Jy/B -> mJy/B


    #positive
    for temp_clump in range(0,Nclump):
        max_peakSN = np.nanmax(max_SNcube[np.where(id_SNcube == temp_clump+1)])
        where_clump = np.array(np.where(max_SNcube == max_peakSN))
        if len(where_clump) > 1:
            len_now = int(len(where_clump[0])/2)
            x_now = where_clump[0,len_now]
            y_now = where_clump[1,len_now]
            z_now = where_clump[2,len_now]
            clfind[temp_clump, 0:3] = np.array([x_now, y_now, z_now])
        else:
            clfind[temp_clump, 0:3] = np.where(max_SNcube == max_peakSN) # position

        clfind[temp_clump, 3] = max_peakSN # max peak SN
        clfind[temp_clump, 4] = data_flux[clfind[temp_clump,0],clfind[temp_clump,1],clfind[temp_clump,2]]
        clfind[temp_clump, 5] = np.size(np.where(id_SNcube == temp_clump+1))/3 # Npixel
        print 'circle(%d %d, 0.5") #text={%d: %.2f sigma, %d sl, %d pix}' % (int(clfind[temp_clump,0]), int(clfind[temp_clump,1]), temp_clump, clfind[temp_clump,3], int(clfind[temp_clump,2]), int(clfind[temp_clump,5]))


    #negative
    for temp_clump in range(0,Nclump_inv):
        max_peakSN = np.nanmin(max_SNcube[np.where(id_SNcube == -1*(temp_clump+1))])
        where_clump = np.array(np.where(max_SNcube == max_peakSN))
        if len(where_clump) > 1:
            len_now = int(len(where_clump[0])/2)
            x_now = where_clump[0,len_now]
            y_now = where_clump[1,len_now]
            z_now = where_clump[2,len_now]
            clfind_inv[temp_clump, 0:3] = np.array([x_now, y_now, z_now])
        else:
            clfind_inv[temp_clump, 0:3] = np.where(max_SNcube == max_peakSN) # position

        clfind_inv[temp_clump, 3] = max_peakSN # max negative peak SN
        clfind_inv[temp_clump, 4] = data_flux[clfind_inv[temp_clump,0],clfind_inv[temp_clump,1],clfind_inv[temp_clump,2]]
        clfind_inv[temp_clump, 5] = np.size(np.where(id_SNcube == -1*(temp_clump+1)))/3 # Npixel
    sorted(clfind, key=itemgetter(3,4,5), reverse=True)
    sorted(clfind_inv, key=itemgetter(3,4))

    del id_SNcube
    print ''


#-------------------------------------------------------------------Step 6 Contamination check
    mystep = 6
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#Step', mystep, step_title[mystep]


    if Nclump > Nclump_inv: tot = Nclump
    else: tot = Nclump_inv

    max_peakSN = np.nanmax(max_SNcube)
    min_peakSN = np.nanmin(max_SNcube)
    del max_SNcube

    if max_peakSN > -1.0*min_peakSN: high = max_peakSN
    else: high = -1.0*min_peakSN

    d_SN = 0.1
    Nbin = int((high - low) / d_SN + 1)

    cum_dif = np.zeros([Nbin,4], dtype = 'int')

    for clump in range(0,Nclump):
        dum = (int)((clfind[clump,3]-low)/d_SN)
        cum_dif[dum,0] += 1
        cum_dif[dum,2] += 1

    for clump in range(0,Nclump_inv):
        dum = (int)((-1.0*clfind_inv[clump,3]-low)/d_SN)
        cum_dif[dum,1] += 1
        cum_dif[dum,3] += 1

    print '#  max SN dist vs clump'

    for SN in range(Nbin-1,1,-1):
        cum_dif[SN-1,0] += cum_dif[SN,0]
        cum_dif[SN-1,1] += cum_dif[SN,1]
        print '%.2f %d %d %d %d' % (SN*d_SN+low, cum_dif[SN,0], cum_dif[SN,1], cum_dif[SN,2], cum_dif[SN,3])


    print '# posi_id, position, max_peakSN, peakflux, contamination'

    for clump in range(0,Nclump):
        dum = (int)((clfind[clump,3]-low)/d_SN)
        clfind[clump, 6] = float(cum_dif[dum,0]/cum_dif[dum,1])
        clfind_res[clump,:,spw_now] = clfind[clump,:]
        print 'spw %d, id %d: (%d, %d, %d) %.2f sigma(peak),  %.2f mJy/B(peak), %.2f parcent contami., %d pix.' \
        % (spw_now, clump, int(clfind[clump,0]), int(clfind[clump,1]), int(clfind[clump,2]), \
        clfind[clump,3], clfind[clump,4], clfind[clump,6]*100, int(clfind[clump,5]))

    print ''


    print '# nega_id, position, max_peakSN, peakflux, contamination'

    for clump in range(0,Nclump_inv):
        clfind_res_inv[clump,:,spw_now] = clfind_inv[clump,:]
        print 'id %d: (%d, %d, %d) %.2f sigma(peak),  %.2f mJy/B(peak), %d pix.' \
        % (clump, int(clfind_inv[clump,0]), int(clfind_inv[clump,1]), int(clfind_inv[clump,2]), \
        clfind_inv[clump,3], clfind_inv[clump,4], int(clfind_inv[clump,5]))
    print ''



#-------------------------------------------------------------------Step 7 List candidates, output results
mystep = 7
if STEP4567 == 0:
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#Step', mystep, step_title[mystep]


    for spw_now in range(num_spw):

        print '# spw, id, position, max_peakSN, peakflux, contamination, Npixel'

        for clump in range(0,NCLUMP_TOT):
            if clfind_res[clump, 5, spw_now] > 0:
                print '%d    %d    %d    %d    %d    %.2f    %.2f    %.2f    %d    ' \
                % (spw_now, clump, int(clfind_res[clump,0,spw_now]), int(clfind_res[clump,1,spw_now]), \
                int(clfind_res[clump,2, spw_now]), clfind_res[clump,3, spw_now], clfind_res[clump,4, spw_now],\
                clfind_res[clump,6, spw_now], int(clfind_res[clump, 5, spw_now]))


        print ''


#-------------------------------------------------------------------Step 8 Completeness check
mystep = 8
if STEP8 == 0:
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print '#Step', mystep, step_title[mystep]


    if not os.path.exists('linesearch/completeness'):
        os.system('rm -rf linesearch/completeness')
    os.system('mkdir linesearch/completeness')

    Set_start = 0
    Set_num =    10

    FWHM_start = 100
    FWHM_end = 1000
    FWHM_bin = 100

    size_start = 5
    size_end = 15
    size_bin = 5
    fac_size = 0.1

    flux_start = 5
    flux_end = 15
    flux_bin = 1
    fac_flux = 0.1

    for spw_now in range(num_spw):
        #seed(Set)
		if STEP123[2] == 0:
			stat = imstat(imagename='linesearch/masked/spw%d.cube.contsub.masked.image' % spw_now,
            axes=[0,1],
            logfile='linesearch/masked/spw%d.cube.contsub.masked.imstat' % spw_now)

		stat_sigma = stat["sigma"]
		chans = chans_array[spw_now]

		start_chan = chans
		end_chan = 0

		for chan_now in xrange(0,chans):
			if  stat_sigma[chan_now] > 0.0:
				if start_chan > chan_now : start_chan = chan_now
				if end_chan < chan_now :   end_chan   = chan_now
		sm_chans = end_chan - start_chan + 1

		V_RESO = np.mean(vel_ori[spw_now])
		sm_min = int(vel_min/V_RESO)
		sm_max = int(vel_max/V_RESO)
		d_sm   = int(Nbin_vel/V_RESO)

		original = imval('linesearch/masked/spw%d.cube.contsub.masked.image' % spw_now,
				chans='%d~%d' % (start_chan, end_chan),
				box='0,0,%d,%d' % (imsize_x-1, imsize_y-1))
		cube = original['data']
		del original

		flux_data = imval('linesearch/original/spw%d.cube.contsub.pb' % spw_now,
                chans='%d~%d' % (start_chan, end_chan),
                box='0,0,%d,%d' % (imsize_x-1, imsize_y-1))

		flux_cube = flux_data['data']
		del flux_data

		for size_now in range(size_start, size_end+size_bin, size_bin):
			size = size_now*fac_size
			source_pixel = size/cell;

			for FWHM in range(FWHM_start, FWHM_end+FWHM_bin, FWHM_bin):
				source_slice = FWHM/V_RESO;

				for flux_now in range(flux_start, flux_end+flux_bin, flux_bin):
					flux = fac_flux*flux_now

					for Set in range(Set_start, Set_start+Set_num):
						cen_i = (int)(imsize_x-1)*np.random.rand()
						cen_j = (int)(imsize_y-1)*np.random.rand()
						cen_k = (int)(start_chan-1 + sm_chans*np.random.rand())

						if cen_i - 0.5*source_pixel < 0:
							minedge_i = 0
						else:
							minedge_i = (int)(cen_i - 0.5*source_pixel)

						if cen_i + 0.5*source_pixel > imsize_x:
							maxedge_i = imsize_x - 1
						else:
							maxedge_i = (int)(cen_i + 0.5*source_pixel)


                        if cen_j - 0.5*source_pixel < 0:
                            minedge_j = 0
                        else:
                            minedge_j = (int)(cen_j - 0.5*source_pixel)

                        if cen_j + 0.5*source_pixel > imsize_y:
                            maxedge_j = imsize_y - 1
                        else:
                            maxedge_j = (int)(cen_j + 0.5*source_pixel)


                        if cen_k - 0.5*source_slice < 0:
                            minedge_k = 0
                        else:
                            minedge_k = (int)(cen_k - 0.5*source_slice)

                        if cen_k + 0.5*source_slice > sm_chans:
                            maxedge_k = sm_chans
                        else:
                            maxedge_k = (int)(cen_k + 0.5*source_slice)


                        print '#        %d %.1f %d %.1f %d %d %d %d %d~%d %d~%d %d~%d' \
                            % (spw_now, size, FWHM, flux, Set, cen_i, cen_j, cen_k,\
                            minedge_i, maxedge_i, minedge_j, maxedge_j, minedge_k, maxedge_k)


                        trim_region = cube[minedge_i:maxedge_i+1,minedge_j:maxedge_j+1,:]
                        smoothed = np.zeros([maxedge_i-minedge_i+1,maxedge_j-minedge_j+1,chans], dtype = 'float')
                        SNcube = np.zeros([maxedge_i-minedge_i+1,maxedge_j-minedge_j+1,chans], dtype = 'float')
                        maxSNcube = np.zeros([maxedge_i-minedge_i+1,maxedge_j-minedge_j+1,chans], dtype = 'float')
                        maxSNmap = np.zeros([maxedge_i-minedge_i+1,maxedge_j-minedge_j+1], dtype = 'float')
                        id_SNcube = np.zeros([maxedge_i-minedge_i+1,maxedge_j-minedge_j+1,chans], dtype = 'int')
                        id_SNmap = np.zeros([maxedge_i-minedge_i+1,maxedge_j-minedge_j+1], dtype = 'int')
                        std = np.zeros([sm_chans], dtype = 'float')

                        for i in range(maxedge_i-minedge_i):
                            for j in range(maxedge_j-minedge_j):
                                for k in range(minedge_k, maxedge_k):
                                    xx = np.power(2.35*(minedge_i+i - cen_i)/source_pixel, 2)#FWHM = 2.35sigma
                                    yy = np.power(2.35*(minedge_j+j - cen_j)/source_pixel, 2)
                                    zz = np.power(2.35*(k - cen_k)/source_slice, 2)
                                    gauss_index= - 0.5 * (xx + yy + zz)
                                    gaussian = flux * np.exp( - gauss_index)
                                    trim_region[i,j,k] = trim_region[i,j,k]\
                                        + gaussian * flux_cube[minedge_i+i,minedge_j+j,k]
                                    #print '%f %d %d %d\n' % (gaussian, minedge_i+i, minedge_j+j, k)

                        for sm_now in range(sm_min,sm_max+d_sm,d_sm):
                            chan_first = int(start_chan +  sm_now/2)
                            chan_last  = int(end_chan - sm_now/2)

                            # smoothing
                            nvalid_chan = end_chan-sm_now-start_chan+1
                            if sm_now == 0: # NO SMOOTHING
                                smoothed[:,:,start_chan:end_chan+1] = cube[:,:,0:sm_chans]

                            else:
                                kernel = np.ones((1,1,sm_now), dtype=np.float64)/float(sm_now)
                                smoothed[:,:,chan_first:nvalid_chan+chan_first] \
                                    = scipy.signal.convolve(trim_region, kernel, mode='same')\
                                        [:,:,sm_now/2:nvalid_chan+sm_now/2]

                            # calc smoothed RMS and normarize
                            std = np.std(smoothed, axis=(0,1))# !--------------------!!!!
                            for chan_center in xrange(chan_first, chan_first+nvalid_chan):
                                SNcube[:,:,chan_center] = smoothed[:,:,chan_center] / std[chan_center]

                        # make max S/N cube
                        if sm_now == sm_min:
                            maxSNcube = SNcube
                        else:
                            maxSNcube = np.fmax(SNcube, maxSNcube)#matching in 3D

                        #max S/N map
                        for chan_now in range(start_chan+1,end_chan+1):#matching in 2D
                            if chan_now == start_chan + 1 :
                                maxSNmap = np.fmax(maxSNcube[:,:,chan_now-1], maxSNcube[:,:,chan_now])
                            else:
                                maxSNmap = np.fmax(maxSNmap, maxSNcube[:,:,chan_now])

                        # check detected or not
                        maxSNcube[(-1.0*low  < maxSNcube) & (maxSNcube < low)] = 0.0
                        maxSNmap[(-1.0*low  < maxSNmap) & (maxSNmap < low)] = 0.0

                        for j in range(maxedge_j-minedge_j):
                            for i in range(maxedge_i-minedge_i):
                                if maxSNmap[i,j] >= low: # when meet new clump at (i+k, j)
                                    if id_SNmap[i-1,j-1] > 0: id_SNmap[i,j] = id_SNmap[i-1,j-1]
                                    elif id_SNmap[i-1,j] > 0: id_SNmap[i,j] = id_SNmap[i-1,j]
                                    elif id_SNmap[i,j-1] > 0: id_SNmap[i,j] = id_SNmap[i,j-1]
                                    elif id_SNmap[i+1,j-1] > 0: id_SNmap[i,j] = id_SNmap[i+1,j-1]
                                    else:
                                        id_SNmap[i, j] = 1

                        temp_id = 0
                        if np.size(np.where(id_SNmap == 1))/2 < NPIXMIN:#!------------------------- /2 ???
                            id_SNmap[id_SNmap == 1] = 0

                        where_source = np.array(np.where(id_SNmap == 1))
                        y_len = len(where_source[0])
                        for k in range(minedge_k, maxedge_k):
                            for y_now in range(0,y_len):
                                i_now = where_source[0,y_now]
                                j_now = where_source[1,y_now]
                                if maxSNcube[i_now,j_now,k] >= low :
                                    if id_SNcube[i_now,j_now,k-1] > 0 :
                                        id_SNcube[i_now,j_now,k] = id_SNcube[i_now,j_now,k-1]
                                    else:
                                        id_SNcube[i_now, j_now,k] = 1

                        if np.size(np.where(id_SNcube == 1))/2 >= NPIXMIN:
                            print 'detected!'


                        del trim_region
                        del smoothed
                        del SNcube
                        del maxSNcube
                        del maxSNmap
                        del id_SNcube
                        del id_SNmap
                        del std

    del cube
    del flux_cube

else:
    print '#Step', mystep, step_title[mystep], ' is skipped.'


print ''
os.system('mkdir linesearch/last')
os.system('mv *.last linesearch/last/')
print 'all steps done.'


import casadef
from macresource import resource_pathname

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



def main ():
    print_and_log('#start read parameter')
    input_params = [
        ### parameters controlling generating datacube
                    'Template/MesurementSet',   #[0] path to the measurement set.
                    1, #[1] epoch. (0 = all epochs, -4 = all epochs - epoch4) (multi epoch is not supported now)
                    2, #[2] Spw (now only spw2 is supported)
                    24, #[3] field 24 (-1 = mosaic) (now only field24 is supported)
                    0, #[4] 1st half or second half of observation (0 = both) 
                    512, #[5] imsize
                    0.12, #[6] source_pixel [arcsec]
                    'J2000 03h32m39.53769 -27d46m43.5673', #[7] phasecenter of the field
                    1, # [8] number of Spws (multi spws is not supported now)
        ### parameters controlling number and property of artificial source
                    0, # [9] start of Set 
                    1, # [10] end of Set
                    100, # [11] FWHM_start in km/s for roop
                    100, # [12] FWHM_end in km/s
                    100, # [13] FWHM_bin in km/s 
                    1,  # [14] size_start in pixel for roop
                    1,  # [15] size_end in pixel
                    1,  # [16] size_bin in pixel (if [14]=[15], pointsource, extend source is not supported.)
                    30, #[17] SN_start [sigma] (peak SN ratio of original freq. resolution) 
                    30, #[18] SN_end [sigma]
                    1, #[19] SN_bin [sigma] (-1 = no for roop)
                    1 #[20] seed of np.rand
                    ] # other parameters will be appended in following steps  
    # a,b,c,d,e = input_params[0], ... , 
    
    check_casa_version()
    generate_template_imagecube(input_params)
    check_parameters(input_params)
    simobserve_simplenoise(input_params)
    imagename_now = input_params[22]+'.flux'  
    ia.open(imagename_now)
    datacube = ia.getchunk()
    ia.close()      

    np.random.seed(seed=input_params[20])
    for FWHM in range(input_params[11], input_params[12]+input_params[13], input_params[13]): 
        for source_pixel in range(input_params[14], input_params[15]+input_params[16], input_params[16]):
            for SN_now in range(input_params[17], input_params[17]+input_params[18], input_params[18]):
                source_slice = FWHM/input_params[35]; #V_RESO
                size = source_pixel*input_params[6] #TBC calc from beam size 
                flux = SN_now*input_params[41] #TBC calc for each slice
                #flux = 0.4
                flux = 0.004
                #flux = 0.0024
                for Set in range(input_params[9], input_params[9]+input_params[10]):
                    cen_i = (int)((input_params[32]-1)*np.random.rand())
                    cen_j = (int)((input_params[33]-1)*np.random.rand())
                    cen_k = (int)((input_params[34] - 1)*np.random.rand())

                    while np.isnan(datacube[cen_i, cen_j, 0, cen_k]):
                        cen_i = (int)((input_params[32]-1)*np.random.rand())
                        cen_j = (int)((input_params[33]-1)*np.random.rand())
                        cen_k = (int)((input_params[34] - 1)*np.random.rand())
                    source_params = [FWHM, size, flux, Set, cen_i, cen_j, cen_k]
                    get_direction_etc(input_params, source_params)
                    
                    simobserve_simplenoise_addsource(input_params, source_params)                 
                    simobserve_addsource_ft_uvsub(input_params, source_params)
                    simobserve_addsource_image(input_params, source_params)            
                    original_addsource_ft_uvsub(input_params, source_params)
                    original_addsource_image(input_params, source_params)            
                                         
                    comparison_data(input_params, source_params)
                    completeness_check(input_params, source_params)    
                    del source_params
    del datacube
    move_logfile()


def print_and_log(message, priority='INFO'):
    print(message)
    casalog.post(message, priority=priority)
    

def check_casa_version():
    version = casadef.casa_version
    print "#You are using " + version

    if (version < '5.1.0'):
        print "#YOUR VERSION OF CASA IS TOO OLD FOR THIS GUIDE."
        print "#PLEASE UPDATE IT BEFORE PROCEEDING."
    else:
        print "#Your version of CASA is appropriate for this guide.\n"


def generate_template_imagecube(input_params):
    print_and_log('\n#generate imagecube.')

    if input_params[1] == 0:         
        vis_now = input_params[0]+'/epochall_hudf_spw%d_scanav_field%d.ms' \
                            % (input_params[2], input_params[3])
        imagename_now ='TemplateImage/epochall_hudf_spw%d_scanav_field%d.clean' \
                        % (input_params[2], input_params[3])
    
    elif input_params[1] > 0:                
        if input_params[4] == 0:
            vis_now = input_params[0]+'/epoch%d_hudf_spw%d_scanav_field%d.ms' \
                            % (input_params[1], input_params[2], input_params[3])
            imagename_now ='TemplateImage/epoch%d_hudf_spw%d_scanav_field%d.clean' \
                            % (input_params[1], input_params[2], input_params[3])
        elif input_params[4] == 1:    
            vis_now = input_params[0]+'/epoch%d_hudf_spw%d_scanav_field%d_1st.ms' \
                            % (input_params[1], input_params[2], input_params[3])
            imagename_now ='TemplateImage/epoch%d_hudf_spw%d_scanav_field%d_1st.clean' \
                            % (input_params[1], input_params[2], input_params[3])
                            
        elif input_params[4] == 2:
            vis_now = input_params[0]+'/epoch%d_hudf_spw%d_scanav_field%d_2nd.ms' \
                            % (input_params[1], input_params[2], input_params[3])
            imagename_now ='TemplateImage/epoch%d_hudf_spw%d_scanav_field%d_2nd.clean' \
                            % (input_params[1], input_params[2], input_params[3])


    if input_params[3] == -1: 
        imagermode_now = 'mosaic'
    else: 
        imagermode_now = 'csclean'

    input_params.append(vis_now) #21
    input_params.append(imagename_now) #22
    input_params.append(imagermode_now) #23

    print_and_log('#vis_now = %s' % input_params[21])
    print_and_log('#imagename_now = %s' % input_params[22])
    print_and_log('#imagermode_now = %s' % input_params[23])

    if not os.path.exists('%s' % input_params[21]):
        raise Exception, "Could not find an input MS: %s"  % input_params[21]


    if os.path.exists('%s.image' % imagename_now):
        print_and_log('#template image cube already exists. skipped\n')
    else:
        if os.path.exists('TemplateImage'):
            os.system('rm -rf TemplateImage')
        os.system('mkdir TemplateImage')
        
        default(clean)    
        clean(vis = vis_now,
              imagename = imagename_now, 
              imagermode=imagermode_now,
              mode='channel',
              weighting='natural',
              niter = 0, 
              imsize= input_params[5],
              cell= '%.2farcsec' % input_params[6],
              phasecenter = input_params[7],
              minpb = 0.01,
              pbcor = True
              )
        
        os.system('mv %s.image %s.pbcor' % (imagename_now, imagename_now))
        immath(
            imagename = ['%s.pbcor' % imagename_now, '%s.flux' % imagename_now],
            expr = 'IM0*IM1',
            outfile = '%s.image' % imagename_now,
            )
        
    
def check_parameters(input_params):
    print_and_log('\n#check parameters.')

    epoch = input_params[1]
    spw = input_params[2]
    field = input_params[3]
    half = input_params[4]
               
    rms_vis = visstat(vis=input_params[21], datacolumn='data', axis='real')['DATA_DESC_ID=0']['rms']

    imagename_now = input_params[22]+'.image'  

    project_1 = 'sim.spw%d.field%d.ep%d.onlysimplenoise' % (spw, field, epoch)
    project_2 = 'sim.spw%d.field%d.ep%d.simplenoise' % (spw, field, epoch)
    project_3 = 'sim.spw%d.field%d.ep%d.ft_uvsub' % (spw, field, epoch)
    project_4 = 'sim.spw%d.field%d.ep%d.fromimage' % (spw, field, epoch)
    project_5 = 'ori.spw%d.field%d.ep%d.ft_uvsub' % (spw, field, epoch)
    project_6 = 'ori.spw%d.field%d.ep%d.fromimage' % (spw, field, epoch)
    
    ptgfile='Template/ptgfiles/sim.ptg%d.txt' % field
    refdate='2014/0%d/%d' % (half, epoch)
    integration='34.272s' #after scan avarage

    if input_params[4] == 0:
        totaltime='2'
        antennalist = "epoch%d_hudf_spw2_scanav_field%d.ms" % (epoch, field)
    elif input_params[4] == 1:    
        totaltime='1' # number without units = time of repeating
        antennalist = "epoch%d_hudf_spw2_scanav_field%d_%dst.ms" % (epoch, field, half)
    elif input_params[4] == 2:
        totaltime='1' # number without units = time of repeating
        antennalist = "epoch%d_hudf_spw2_scanav_field%d_%dnd.ms" % (epoch, field, half)

    ia.open(imagename_now)
    datacube = ia.getchunk()
    imsize_x = ia.shape()[0]
    imsize_y = ia.shape()[1]
    chans = ia.shape()[3]   
    V_RESO = 20.52     #TBC: get value from data km/s
    inwidth = 15.625      #TBC: get value from data [MHz]
    beam1 =  0.467       #TBC: get value from data #arcsec
    beam2 =  0.360       #TBC: get value from data
    incenter = 228.1776275635 #TBC: get value from data [GHz]
    ia.close()      

    std_ori = np.nanstd(datacube, axis=(0,1,2))     
    fac_flux = np.nanmedian(std_ori)  #mJy/beam
    del datacube

    input_params.append(rms_vis) #24 
    input_params.append(project_1) #25
    input_params.append(project_2) #26
    input_params.append(ptgfile)  #27
    input_params.append(refdate)  #28
    input_params.append(integration)  #29
    input_params.append(totaltime)  #30
    input_params.append(antennalist)  #31
    input_params.append(imsize_x) #32
    input_params.append(imsize_y) #33
    input_params.append(chans) #34
    input_params.append(V_RESO) #35
    input_params.append(inwidth) #36
    input_params.append(beam1) #37
    input_params.append(beam2) #38
    input_params.append(incenter) #39
    input_params.append(std_ori) #40
    input_params.append(fac_flux) #41

    input_params.append(project_3) #42
    input_params.append(project_4) #43
    input_params.append(project_5) #44
    input_params.append(project_6) #45

    print_and_log('#rms_vis(ori) = %f Jy' % input_params[24]) # check unit
    print_and_log('#rms_image(ori) = %f Jy/B' %  input_params[41])


def simobserve_simplenoise(input_params):
    print_and_log('\n#simobserve without source.')

    if os.path.exists('%s' % input_params[25]):
        print_and_log('#directory %s already exists. skipped.' % input_params[25])
    else:

        newim = ia.newimagefromimage( infile='%s.flux' % input_params[22], outfile='TemplateImage/zero.image',overwrite=True)        
        datacube = newim.getchunk()
        datacube[~np.isnan(datacube)] =0.0
        newim.putchunk(datacube)
        newim.close()
        del datacube
        
        os.system('rm -rf TemplateImage/zero.fits')   
        exportfits(imagename='TemplateImage/zero.image', fitsimage='TemplateImage/zero.fits')
        os.system('rm -rf TemplateImage/zero.image')   
        
        if input_params[1] == 0:
            vis_now = []
            
            for epoch in range(1,14):
                if epoch != 4:
                    vis_temp = input_params[0]+'/epoch%d_hudf_spw%d_scanav_field%d.ms' \
                            % (epoch, input_params[2], input_params[3])
                    rms = visstat(vis=vis_temp, datacolumn='data', axis='real')['DATA_DESC_ID=0']['rms']
                    refdate_now='2014/01/%d' % epoch
                    antennalist_now = "epoch%d_hudf_spw2_scanav_field%d.ms" % (epoch, input_params[3])
                    simobserve(
                        skymodel = 'TemplateImage/zero.fits',
                        inbright = '0.000Jy/pixel', 
                        indirection= '%s' % input_params[7],
                        incell = '%.2farcsec' % input_params[6],
                        incenter = '%fGHz' % input_params[39],
                        inwidth = '%fMHz' % input_params[36], 
                        project = '%s' % input_params[25], 
                        #
                        setpointings=False,
                        ptgfile='%s' % input_params[27], 
                        #
                        obsmode='int',
                        refdate='%s' % refdate_now, 
                        integration='%s' % input_params[29], 
                        direction= '%s'% input_params[7], # center of field 24
                        totaltime='2', # number without units = time of repeating
                        antennalist = "Template/configurationfiles/%s.cfg" %  antennalist_now, 
                        overwrite = True
                        )
                    
                    vis_temp = '%s/%s.%s.ms' % (input_params[25], input_params[25], antennalist_now)
                    vis_now.append(vis_temp)
                    sm.openfromms(vis_temp)
                    sm.setnoise(mode='simplenoise', simplenoise="%fJy" % rms)
                    sm.corrupt()
                    sm.done()               
                      
            concat(vis = vis_now, concatvis='%s/%s.ms' % (input_params[25], input_params[25]))
            del vis_now 
            vis_now = '%s/%s.ms' % (input_params[25], input_params[25])
            

        else:
            simobserve(
                skymodel = 'TemplateImage/zero.fits',
                inbright = '0.000Jy/pixel', 
                indirection= '%s' % input_params[7],
                incell = '%.2farcsec' % input_params[6],
                incenter = '%fGHz' % input_params[39],
                inwidth = '%fMHz' % input_params[36], 
                project = '%s' % input_params[25], 
                #
                setpointings=False,
                ptgfile='%s' % input_params[27], 
                #
                obsmode='int',
                refdate='%s' % input_params[28], 
                integration='%s'% input_params[29], #after scan avarage
                direction= '%s'% input_params[7], # center of field 24
                totaltime='%s'% input_params[30], # number without units = time of repeating
                antennalist = "Template/configurationfiles/%s.cfg" % input_params[31], 
                overwrite = True
                )
                    
            vis_now = '%s/%s.%s.ms' % (input_params[25], input_params[25], input_params[31])
            sm.openfromms(vis_now)
            sm.setnoise(mode='simplenoise', simplenoise="%fJy" % input_params[24])
            sm.corrupt()
            sm.done()


        imagename_now = '%s/%s' % (input_params[25], input_params[25])
        #print imagename_now

        default(clean)    
        clean(vis = vis_now,
              imagename = imagename_now, 
              imagermode=input_params[23],
              mode='channel',
              weighting='natural',
              niter = 0, 
              imsize= input_params[5],
              cell= '%.2farcsec' % input_params[6],
              phasecenter = input_params[7],
              minpb = 0.01,
              pbcor = True
              )
        
        os.system('mv %s.image %s.pbcor' % (imagename_now, imagename_now))
        immath(
            imagename = ['%s.pbcor' % imagename_now, '%s.flux' % imagename_now],
            expr = 'IM0*IM1',
            outfile = '%s.image' % imagename_now,
            )

        ia.open('%s.image' % imagename_now)
        datacube = ia.getchunk()
        ia.close()      

        std_ori = np.nanstd(datacube, axis=(0,1,2))     
        fac_flux = np.nanmedian(std_ori)  #mJy/beam
        del datacube
        rms = visstat(vis=vis_now, datacolumn='data', axis='real')['DATA_DESC_ID=0']['rms']

        print_and_log('#rms_vis(sim) = %f Jy' % rms) # check unit
        print_and_log('#rms_image(sim) = %f Jy/B' %  fac_flux)


def get_direction_etc(input_params, source_params):
    print_and_log('\n#get direction etc.')

    #TBC convert pixel to RADEC
    source_params[4] = 256
    source_params[5] = 281 #3arcmin offset
    #source_params[5] = 306 #6
    #source_params[5] = 331 #9 
    source_params[6] = 60
    indirection= 'J2000 03h32m39.465 -27d46m40.687', # source position  
    #indirection= 'J2000 03h32m39.465 -27d46m37.687', # source position  
    #indirection= 'J2000 03h32m39.465 -27d46m34.687', # source position  

    
    imagename_now = input_params[22]+'.flux'  
    ia.open(imagename_now)
    datacube = ia.getchunk()
    ia.close()      
    
    source_pixel = (int)(source_params[1]/input_params[6])
    source_slice = (int)(source_params[0]/input_params[35])    
    if source_params[4] - source_pixel < 0:
        minedge_i = 0
    else:
        minedge_i = (int)(source_params[4] - source_pixel) 
                                       
    if source_params[4] + source_pixel > input_params[32]-1:
        maxedge_i = input_params[32] - 1
    else:
        maxedge_i = (int)(source_params[4] + source_pixel)
    
    if source_params[5] - source_pixel < 0:
        minedge_j = 0
    else:
        minedge_j = (int)(source_params[5] - source_pixel)
            
    if source_params[5] + source_pixel > input_params[33]-1:
        maxedge_j = input_params[33] - 1
    else:
        maxedge_j = (int)(source_params[5] + source_pixel)

    if source_params[6] - source_slice < 0:
        minedge_k = 0
    else:
        minedge_k = (int)(source_params[6] - source_slice)
            
    if source_params[6] + source_slice > input_params[34]-1:
        maxedge_k = input_params[34]-1
    else:
        maxedge_k = (int)(source_params[6] + source_slice)

    source_params.append(indirection) #7
    source_params.append(minedge_i) #8
    source_params.append(maxedge_i) #9
    source_params.append(minedge_j) #10
    source_params.append(maxedge_j) #11
    source_params.append(minedge_k) #12
    source_params.append(maxedge_k) #13
        
    #make model array 
    print_and_log('\n#make model data.')

    os.system('rm -rf TemplateImage/line.image')
    newim = ia.newimagefromimage( infile='%s.image' % input_params[22], 
                                  outfile='TemplateImage/line.image',
                                  overwrite=True)        
    datacube = newim.getchunk()
    datacube[~np.isnan(datacube)] =0.0
    
    flux = source_params[2]
    if source_pixel == 1:
        for k in range(maxedge_k-minedge_k+1):
            if not np.isnan(datacube[source_params[4], source_params[5], 0, k+minedge_k]):
                xx = 0.0
                yy = 0.0
                zz = np.power(2.35*(minedge_k+k - source_params[6])/source_slice, 2)                            
                gauss_index= - 0.5 * (xx + yy + zz)
                gaussian = flux * np.exp(gauss_index)
                datacube[source_params[4],source_params[5],:,k+minedge_k] += gaussian
                #print ('%f %f %d %d %d' % (gauss_index, gaussian, source_params[4], source_params[5], minedge_k+k)) 
    else:
        for k in range(maxedge_k-minedge_k+1):
            for i in range(maxedge_i-minedge_i+1):
                for j in range(maxedge_j-minedge_j+1):
                    if not np.isnan(datacube[i+minedge_i, j+minedge_j, 0, k+minedge_k]):
                        xx = np.power(2.35*(minedge_i+i - source_params[4])/source_pixel, 2)#FWHM = 2.35sigma
                        yy = np.power(2.35*(minedge_j+j - source_params[5])/source_pixel, 2)
                        zz = np.power(2.35*(minedge_k+k - source_params[6])/source_slice, 2)                            
                        gauss_index= - 0.5 * (xx + yy + zz)
                        gaussian = flux * np.exp(gauss_index)
                        datacube[i+minedge_i,j+minedge_j,:,k+minedge_k] += gaussian
                        #print ('%f %f %d %d %d' % (gauss_index, gaussian, minedge_i+i, minedge_j+j, minedge_k+k)) 

    newim.putchunk(datacube)
    newim.close()
    del datacube
    os.system('rm -rf TemplateImage/line.fits')
    exportfits(imagename='TemplateImage/line.image', fitsimage='TemplateImage/line.fits')

    #convolution 
    print_and_log('\n#convolution.')
    if not os.path.exists('%s_conv.PBresponsed' % input_params[22]):
        ia.open('%s.flux' % input_params[22])
        flux_cube = ia.getchunk()
        ia.close()        
        '''
        myhead   =  imhead('%s.image' % input_params[22], mode='list')
        for k in range(chans)
            bmaj     = myhead['perplanebeams']['*'+str(k)]['major']['value']
            bmin     = myhead['perplanebeams']['*'+str(k)]['minor']['value']
            bpa      = myhead['perplanebeams']['*'+str(k)]['positionangle']['value']
            # how can I put beamsize for each channels?
        ia.restoringbeam(). 
        imsmooth 
        '''
        ia.open('TemplateImage/line.image')
        ia.convolve2d(outfile  = '%s_conv.intrinsic' % input_params[22], 
                    axes=[0,1], type='gauss',  major='0.489arcsec', minor='0.449arcsec', pa='-56.961deg')
        ia.done()
        newim = ia.newimagefromimage( infile='%s_conv.intrinsic' % input_params[22], 
                                      outfile='%s_conv.PBresponsed' % input_params[22], overwrite=True)
        line_array = newim.getchunk()        
        line_array = line_array * flux_cube
        newim.putchunk(line_array)
        newim.close()
        del line_array, flux_cube

    else:
        print_and_log('convolved model file already exists. skipped.')
        
    print ('#source_param: %dkm/s, %farcsec, %fJy/B, %dSet, [%d %d %d]' % (source_params[0],source_params[1],source_params[2],source_params[3],source_params[4],source_params[5],source_params[6]))
    print ('%d %d' % (source_pixel, source_slice))
    print ('%d %d %d %d %d %d' % (minedge_i, maxedge_i, minedge_j, maxedge_j, minedge_k, maxedge_k))
    
    
def simobserve_simplenoise_addsource(input_params, source_params):
    print_and_log('\n#simobserve with source.')
   
    if os.path.exists('%s' % input_params[26]):
        print_and_log('#directory %s already exists. skipped.' % input_params[26])
    else:

        if input_params[1] == 0:
            vis_now = []
            vis_now_zero = []
            
            for epoch in range(1,14):
                if epoch != 4:
                    vis_temp = input_params[0]+'/epoch%d_hudf_spw%d_scanav_field%d.ms' \
                            % (epoch, input_params[2], input_params[3])
                    rms = visstat(vis=vis_temp, datacolumn='data', axis='real')['DATA_DESC_ID=0']['rms']
                    refdate_now='2014/01/%d' % epoch
                    antennalist_now = "epoch%d_hudf_spw2_scanav_field%d.ms" % (epoch, input_params[3])
                    simobserve(
                        skymodel = 'TemplateImage/line.fits',
                        inbright = '%fJy/pixel' % source_params[2], 
                        indirection= '%s' % source_params[7],
                        incell = '%.2farcsec' % input_params[6],
                        incenter = '%fGHz' % input_params[39],
                        inwidth = '%fMHz' % input_params[36], 
                        project = '%s' % input_params[26], 
                        #
                        setpointings=False,
                        ptgfile='%s' % input_params[27], 
                        #
                        obsmode='int',
                        refdate='%s' % refdate_now, 
                        integration='%s' % input_params[29], 
                        direction= '%s'% input_params[7], # center of field 24
                        totaltime='2', # number without units = time of repeating
                        antennalist = "Template/configurationfiles/%s.cfg" %  antennalist_now, 
                        overwrite = True
                        )
                    
                    vis_temp = '%s/%s.%s' % (input_params[25], input_params[25], antennalist_now)
                    os.system('cp -r %s.ms %s' % (vis_temp, vis_temp))  
                    sm.openfromms(vis_temp)
                    sm.setnoise(mode='simplenoise', simplenoise="%fJy" % rms)
                    sm.corrupt()
                    sm.done()               
                    vis_now.append(vis_temp)
                    vis_temp ='%s/%s.%s.ms' % (input_params[25], input_params[25], antennalist_now)
                    vis_now_zero.append(vis_temp)
                      
            concat(vis = vis_now, concatvis='%s/%s.ms' % (input_params[25], input_params[25]))
            concat(vis = vis_now_zero, concatvis='%s/%s.ms.ms' % (input_params[25], input_params[25]))

            del vis_now, vis_now_zero
            vis_now = '%s/%s.ms' % (input_params[25], input_params[25])

        else:

            simobserve(
                skymodel = 'TemplateImage/line.fits',
                inbright = '%fJy/pixel' % source_params[2], 
                indirection= '%s' % source_params[7],
                incell = '%.2farcsec' % input_params[6],
                incenter = '%fGHz' % input_params[39],
                inwidth = '%fMHz' % input_params[36], 
                project = '%s' % input_params[26], 
                #
                setpointings=False,
                ptgfile='%s' % input_params[27], 
                #
                obsmode='int',
                refdate='%s' % input_params[28], 
                integration='%s'% input_params[29], #after scan avarage
                direction= '%s'% input_params[7], # center of field 24
                totaltime='%s'% input_params[30], # number without units = time of repeating
                antennalist = "Template/configurationfiles/%s.cfg" % input_params[31], 
                overwrite = True
                )

            vis_now = '%s/%s.%s' % (input_params[26], input_params[26], input_params[31])
            os.system('cp -r %s.ms %s' % (vis_now, vis_now))  

            sm.openfromms(vis_now)
            sm.setnoise(mode='simplenoise', simplenoise="%fJy" % input_params[24])
            sm.corrupt()
            sm.done()

        imagename_now = '%s/%s_zeronoise' % (input_params[26], input_params[26])
        default(clean)    
        clean(vis = '%s.ms' % vis_now,
              imagename = imagename_now, 
              imagermode=input_params[23],
              mode='channel',
              weighting='natural',
              niter = 0, 
              imsize= input_params[5],
              cell= '%.2farcsec' % input_params[6],
              phasecenter = input_params[7],
              minpb = 0.01,
              pbcor = True
              )
        
        os.system('mv %s.image %s.pbcor' % (imagename_now, imagename_now))
        immath(
            imagename = ['%s.pbcor' % imagename_now, '%s.flux' % imagename_now],
            expr = 'IM0*IM1',
            outfile = '%s.image' % imagename_now,
            )
        

        imagename_now = '%s/%s' % (input_params[26], input_params[26])
        default(clean)    
        clean(vis = vis_now,
              imagename = imagename_now, 
              imagermode=input_params[23],
              mode='channel',
              weighting='natural',
              niter = 0, 
              imsize= input_params[5],
              cell= '%.2farcsec' % input_params[6],
              phasecenter = input_params[7],
              minpb = 0.01,
              pbcor = True
              )
        
        os.system('mv %s.image %s.pbcor' % (imagename_now, imagename_now))
        immath(
            imagename = ['%s.pbcor' % imagename_now, '%s.flux' % imagename_now],
            expr = 'IM0*IM1',
            outfile = '%s.image' % imagename_now,
            )
        
        ia.open('%s.image' % imagename_now)
        datacube = ia.getchunk()
        ia.close()      

        std_ori = np.nanstd(datacube, axis=(0,1,2))     
        fac_flux = np.nanmedian(std_ori)  #mJy/beam
        del datacube
        rms = visstat(vis=vis_now, datacolumn='data', axis='real')['DATA_DESC_ID=0']['rms']
        print_and_log('#rms_vis(sim+source) = %f Jy' % rms) # check unit
        print_and_log('#rms_image(sim+source) = %f Jy/B' %  fac_flux)
          
    
    
def simobserve_addsource_ft_uvsub(input_params, source_params):

    print_and_log('\n#simobserve ft and uvsub.')
    
    if os.path.exists('%s' % input_params[42]):
        print_and_log('#directory %s already exists. skipped.\n' % input_params[42])
    else:
        os.system('mkdir %s' % input_params[42])

        if input_params[1] == 0:         
            vis_prev = '%s/%s.ms' % (input_params[25], input_params[25])
        else:
            vis_prev = '%s/%s.%s.ms' % (input_params[25], input_params[25], input_params[31])
        
        vis_now = '%s/%s.%s.ms' % (input_params[42], input_params[25], input_params[31])
        os.system('cp -r %s %s' % (vis_prev, input_params[42]))
        
        source_pixel = (int)(source_params[1]/input_params[6])
        source_slice = (int)(source_params[0]/input_params[35])
        flux = source_params[2]
        if source_pixel == 1: shape_now = 'point'
        else: shape_now = 'gauss'        

        imagename_now = input_params[22]+'.flux'  
        ia.open(imagename_now)
        datacube = ia.getchunk()
        ia.close()      
        i = source_params[4]
        j = source_params[5] 

        cl.addcomponent(
            dir='%s' % source_params[7], 
            flux=0.000,
            fluxunit='Jy', 
            shape=shape_now
            ) 
        cl.rename('%s/model00.cl' % input_params[42])
        cl.done()  
                     
        ft(vis=vis_now, 
            spw='0:0~%d' % (input_params[12]),
            complist='%s/model00.cl' % (input_params[42]),
            usescratch=True)            
        ft(vis=vis_now, 
            spw='0:%d~%d' % (input_params[13]+1,input_params[34]-1),
            complist='%s/model00.cl' % (input_params[42]),
            usescratch=True)            
        
        for k_now in range(source_params[12],source_params[13]+1): 
            if not np.isnan(datacube[i, j, 0, k_now]):   
                zz = np.power(2.35*(k_now - source_params[6])/source_slice, 2)
                gauss_index= - 0.5 * zz
                gaussian = flux * np.exp(gauss_index)    
                if source_pixel == 1:
                    cl.addcomponent(
                        dir='%s' % source_params[7], 
                        flux=gaussian,
                        fluxunit='Jy', 
                        shape=shape_now
                        #,majoraxis=str(major_ax)+'arcsec', #TBA, for extended source
                        #minoraxis=str(minor_ax)+'arcsec',  #TBA, for extended source
                        #positionangle=str(angle)+'deg')    #TBA, for extended source
                        ) 

                    cl.rename('%s/model%d.cl' % (input_params[42], k_now))
                    cl.done()  
                    
                    ft(vis=vis_now, 
                       spw='0:%d~%d' % (k_now, k_now),
                       complist='%s/model%d.cl' % (input_params[42], k_now),
                       usescratch=True)            
        del datacube
 
        uvsub(vis=vis_now,reverse=True)
            
        imagename_now = '%s/%s_ft' % (input_params[42], input_params[25])
        default(clean)    
        clean(vis = vis_now,
              imagename = imagename_now, 
              imagermode=input_params[23],
              mode='channel',
              weighting='natural',
              niter = 0, 
              imsize= input_params[5],
              cell= '%.2farcsec' % input_params[6],
              phasecenter = input_params[7],
              minpb = 0.01,
              pbcor = True
              )
        
        os.system('mv %s.image %s.pbcor' % (imagename_now, imagename_now))
        immath(
            imagename = ['%s.pbcor' % imagename_now, '%s.flux' % imagename_now],
            expr = 'IM0*IM1',
            outfile = '%s.image' % imagename_now,
            )


def simobserve_addsource_image(input_params, source_params):
    print_and_log('\n#simobserve image.')
    if os.path.exists('%s' % input_params[43]):
        print_and_log('#directory %s already exists. skipped.\n' % input_params[43])
    else:
        imagename_now = '%s/%s_newim.image' % (input_params[43], input_params[31] )
        os.system('mkdir %s' % input_params[43])

        ia.open('%s_conv.PBresponsed' % input_params[22])
        line_array = ia.getchunk()
        ia.close()
        
        newim = ia.newimagefromimage( infile='%s/%s.image' % (input_params[25], input_params[25]), 
                                      outfile=imagename_now,overwrite=True)
        datacube = newim.getchunk()
        datacube += line_array  
        newim.putchunk(datacube)
        newim.close()
        del datacube, line_array        
        
        
def original_addsource_ft_uvsub(input_params, source_params):
    print_and_log('\n#original ft and uvsub.')
    if os.path.exists('%s' % input_params[44]):
        print_and_log('#directory %s already exists. skipped.\n' % input_params[44])
    else:
        os.system('mkdir %s' % input_params[44])
        
        vis_prev = '%s' % (input_params[21])
        os.system('cp -r %s %s' % (vis_prev, input_params[44]))
        if input_params[1] > 0:                
            if input_params[4] == 0:
                vis_now = input_params[44]+'/epoch%d_hudf_spw%d_scanav_field%d.ms' \
                    % (input_params[1], input_params[2], input_params[3])
            elif input_params[4] == 1:    
                vis_now = input_params[44]+'/epoch%d_hudf_spw%d_scanav_field%d_1st.ms' \
                            % (input_params[1], input_params[2], input_params[3])
            elif input_params[4] == 2:
                vis_now = input_params[44]+'/epoch%d_hudf_spw%d_scanav_field%d_2nd.ms' \
                            % (input_params[1], input_params[2], input_params[3])
        
        source_pixel = (int)(source_params[1]/input_params[6])
        source_slice = (int)(source_params[0]/input_params[35])    
        flux = source_params[2]
        if source_pixel == 1: shape_now = 'point'
        else: shape_now = 'gauss'        

        imagename_now = input_params[22]+'.flux'  
        ia.open(imagename_now)
        datacube = ia.getchunk()
        ia.close()      
        i = source_params[4]
        j = source_params[5]

        cl.addcomponent(
            dir='%s' % source_params[7], 
            flux=0.000,
            fluxunit='Jy', 
            shape=shape_now
            ) 
        cl.rename('%s/model00.cl' % input_params[44])
        cl.done()  
        
        ft(vis=vis_now, 
            spw='0:0~%d' % (input_params[12]),
            complist='%s/model00.cl' % (input_params[44]),
            usescratch=True)            
        ft(vis=vis_now, 
            spw='0:%d~%d' % (input_params[13]+1,input_params[34]-1),
            complist='%s/model00.cl' % (input_params[44]),
            usescratch=True)            
         
        for k_now in range(source_params[12],source_params[13]+1): 
            if not np.isnan(datacube[i, j, 0, k_now]):   
                zz = np.power(2.35*(k_now - source_params[6])/source_slice, 2)
                gauss_index= - 0.5 * zz
                gaussian = flux * np.exp(gauss_index)    
                if source_pixel == 1:
                    cl.addcomponent(
                        dir='%s' % source_params[7], 
                        flux=gaussian,
                        fluxunit='Jy', 
                        shape=shape_now
                        #,majoraxis=str(major_ax)+'arcsec', #TBA, for extended source
                        #minoraxis=str(minor_ax)+'arcsec',  #TBA, for extended source
                        #positionangle=str(angle)+'deg')    #TBA, for extended source
                        ) 

                    cl.rename('%s/model%d.cl' % (input_params[44], k_now))
                    cl.done()  
                    
                    ft(vis=vis_now, 
                       spw='0:%d~%d' % (k_now, k_now),
                       complist='%s/model%d.cl' % (input_params[44], k_now),
                       usescratch=True)            
                        
        del datacube
 
        uvsub(vis=vis_now,reverse=True)

        imagename_now = '%s/%s_ft' % (input_params[44], input_params[31])

        default(clean)    
        clean(vis = vis_now,
              imagename = imagename_now, 
              imagermode=input_params[23],
              mode='channel',
              weighting='natural',
              niter = 0, 
              imsize= input_params[5],
              cell= '%.2farcsec' % input_params[6],
              phasecenter = input_params[7],
              minpb = 0.01,
              pbcor = True
              )
        
        os.system('mv %s.image %s.pbcor' % (imagename_now, imagename_now))
        immath(
            imagename = ['%s.pbcor' % imagename_now, '%s.flux' % imagename_now],
            expr = 'IM0*IM1',
            outfile = '%s.image' % imagename_now,
            )


def original_addsource_image(input_params, source_params):
    print_and_log('\n#original image.')

    if os.path.exists('%s' % input_params[45]):
        print_and_log('#directory %s already exists. skipped.\n' % input_params[45])
    else:
        imagename_now = '%s/%s_newim.image' % (input_params[45], input_params[31] )
        os.system('mkdir %s' % input_params[45])

        ia.open('%s_conv.PBresponsed' % input_params[22])
        line_array = ia.getchunk()
        ia.close()
                    
        newim = ia.newimagefromimage( infile='%s.image' % input_params[22], 
                                      outfile=imagename_now,overwrite=True)
        datacube = newim.getchunk()
        datacube += line_array  
        newim.putchunk(datacube)
        newim.close()
        del datacube, line_array        


def comparison_data(input_params, source_params):
    print 'TBA'
    '''
    if os.path.exists('Comparison'):
        print_and_log('#directory Comparison already exists. skipped.\n')
    else:    
        trim_region = datacube[minedge_i:maxedge_i+1,minedge_j:maxedge_j+1,:,:]
        
        ia.fromarray(pixels=trim_region,
        outfile='newimagefromimage/spw%d.trim.%.1f.%d.%.1f.%d.image' % (spw_now, flux*1e3, FWHM, size, Set))
        ia.done()
        # uvmodelfit  and subtraction 
        # calc residual 
    '''    
        

        
def completeness_check(input_params, source_params):
    print 'TBA'   
    #put your source finding method here! #

def move_logfile():
    if not os.path.exists('last'):
        os.system('mkdir last')
    os.system('mv *.last last/')

if __name__ == '__main__':
    main()

# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

from __future__ import division
from astropy.io import fits
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt

import requests
import random
import datetime
import os,urllib

gzpath = '/Users/willettk/Astronomy/Research/GalaxyZoo'
min_pixelscale = 0.10

def get_nsa_images(highz=True):

    # Should be 145,155 galaxies (as of 26 May 2015)

    if highz:
        version = '1_0_0'
    else:
        version = '0_1_2'

    nsa = fits.getdata('%s/decals/fits/nsa_v%s.fits' % (gzpath,version),1)

    return nsa

def get_decals_bricks(dr='1'):

    # Return the catalog of DECaLS image bricks

    '''
    Include only bricks with images in r, g, and z-bands.
    No constraints on catalog entry for DR2, oddly.
    '''

    bricks_all = fits.getdata('%s/decals/fits/decals-bricks-dr%s.fits' % (gzpath,dr),1)

    if dr == '1':
        has_g = bricks_all['has_image_g']
        has_r = bricks_all['has_image_r']
        has_z = bricks_all['has_image_z']
        has_catalog = bricks_all['has_catalog']

        bricks = bricks_all[has_g & has_r & has_z & has_catalog]

    elif dr == '2':

        has_g = bricks_all['nobs_max_g'] > 0
        has_r = bricks_all['nobs_max_r'] > 0
        has_z = bricks_all['nobs_max_z'] > 0

        bricks = bricks_all[has_g & has_r & has_z]

    return bricks

def find_matching_brick(gal,bricks):

    # Bricks are roughly 0.25 deg x 0.25 deg on each side
    # Three sets of coordinates:
    # ra,dec = center of brick
    # ra1,dec1 = lower right corner of brick
    # ra2,dec2 = upper left corner of brick

    # Find only bricks that could match in RA,dec
    ragal,decgal = gal['RA'],gal['DEC']

    radiff_max = 0.30
    decdiff_max = 0.25

    ramatch = (bricks['ra1'] < ragal) & (bricks['ra2'] >= ragal)
    decmatch = (bricks['dec1'] < decgal) & (bricks['dec2'] >= decgal)

    coomatch = (ramatch & decmatch)
    nmatch = sum(coomatch)

    # Should only be 1 brick per galaxy max, but check
    '''
    if nmatch > 1:
        print 'Attention: %i matches for NSA %s (%.3f,%.3f)' % (nmatch,gal['IAUNAME'],ragal,decgal)
        for idx,brick in enumerate(bricks[coomatch]):
            print '\tBrick #%i: RA from %.3f to %.3f, dec from %.3f to %.3f' % (idx,brick['ra1'],brick['ra2'],brick['dec1'],brick['dec2'])
    '''

    return nmatch,coomatch

def run_all_bricks(nsa,bricks,dr,run_to=-1):

    # Create a matched catalogue of all NSA sources that have grz imaging in DECaLS.
    # Fed into run_NSA() to download FITS files and actually make images.

    '''
    Make this routine somewhat quicker by first eliminating everything in the NSA catalog
    outside the observed RA/dec range of the DECaLS bricks.
    '''

    '''
    Rough limits of DR1 in RA
    Three blocks:
      7h to 11h
      21h to 1h
      15h to 18h
    '''

    brick_maxdec = max(bricks['dec2'])
    brick_mindec = min(bricks['dec1'])

    #ralim = ((nsa['RA'] > 15/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 1/24. * 360) | (nsa['RA'] > 21/24. * 360) | (nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 11/24. * 360)
    #declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

    # Rough limits of DR2 in RA

    ralim = ((nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 3/24. * 360) | (nsa['RA'] > 21/24. * 360) 
    declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

    total_matches = 0
    multi_matches = 0
    decals_indices = np.zeros(len(nsa),dtype=bool)
    bricks_indices = []

    for idx,gal in enumerate(nsa[:run_to]):
        if (declim & ralim)[idx]:
            nm,coomatch = find_matching_brick(gal,bricks)
            if nm > 0:
                total_matches += 1
                decals_indices[idx] = True
                bricks_indices.append(coomatch.argmax())
            if nm > 1:
                multi_matches += 1
        if idx % 100 == 0 and idx > 0:
            print '%i galaxies searched, %i matches so far' % (idx,total_matches)

    print '%i total matches between NASA-Sloan Atlas and DECaLS DR1' % total_matches
    print '%i galaxies had matches in more than one brick' % multi_matches

    nsa_table = Table(nsa)
    nsa_decals = nsa_table[decals_indices]

    assert len(nsa_decals) == len(bricks_indices), \
        "Length of nsa_decals (%i) and bricks_indices (%i) must match" % (len(nsa_decals),len(bricks_indices))

    for bc in bricks.columns:
        nsa_decals[bc.name] = np.array(len(nsa_decals),dtype=bricks[bc.name].dtype)
        for ii in range(len(nsa_decals)):
            nsa_decals[ii][bc.name] = bricks[bricks_indices[ii]][bc.name]

    # Write to file
    # Check what version of the NSA is being used and set string variable below
    nsa_version = '1_0_0'

    outfilename = '%s/decals/fits/nsa_v%s_decals_dr%s.fits' % (gzpath,nsa_version,dr)
    nsa_decals.write(outfilename,overwrite=True)

    return nsa_decals

def get_skyserver_fits(gal,fitspath='%s/decals/fits/nsa' % gzpath,dr='1',remove_multi_fits=True):

    # Get FITS

    galname = gal['IAUNAME']
    params = urllib.urlencode({'ra':gal['RA'],'dec':gal['DEC'],'pixscale':max(min(gal['PETROTH50']*0.04,gal['PETROTH90']*0.02),min_pixelscale),'size':424})
    if dr == '1':
        url = "http://imagine.legacysurvey.org/fits-cutout-decals-dr1?%s" % (params)
    elif dr == '2':
        url = "http://legacysurvey.org/viewer/fits-cutout-decals-dr2?%s" % (params)
    urllib.urlretrieve(url, "%s/%s.fits" % (fitspath,galname))

    # Write multi-plane FITS images to separate files for each band

    data,hdr = fits.getdata("%s/%s.fits" % (fitspath,galname),0,header=True)
    for idx,plane in enumerate('grz'):
        hdr_copy = hdr.copy()
        hdr_copy['NAXIS'] = 2
        hdr_copy['FILTER'] = '%s       ' % plane
        for badfield in ('BANDS','BAND0','BAND1','BAND2','NAXIS3'):
            hdr_copy.remove(badfield)
        fits.writeto("%s/%s_%s.fits" % (fitspath,galname,plane),data[idx,:,:],hdr_copy,clobber=True)


    if remove_multi_fits:
        os.remove("%s/%s.fits" % (fitspath,galname))
    del data,hdr

    return None

def dstn_rgb(imgs, bands, mnmx=None, arcsinh=None, scales=None, desaturate=False):

    '''
    Given a list of images in the given bands, returns a scaled RGB
    image.

    *imgs*  a list of numpy arrays, all the same size, in nanomaggies
    *bands* a list of strings, eg, ['g','r','z']
    *mnmx*  = (min,max), values that will become black/white *after* scaling.
           Default is (-3,10)
    *arcsinh* use nonlinear scaling as in SDSS
    *scales*

    Returns a (H,W,3) numpy array with values between 0 and 1.
    '''
    bands = ''.join(bands)

    grzscales = dict(g = (2, 0.0066),
                      r = (1, 0.01385),
                      z = (0, 0.025),
                      )

    if scales is None:
        if bands == 'grz':
            scales = grzscales
        elif bands == 'urz':
            scales = dict(u = (2, 0.0066),
                          r = (1, 0.01),
                          z = (0, 0.025),
                          )
        elif bands == 'gri':
            scales = dict(g = (2, 0.002),
                          r = (1, 0.004),
                          i = (0, 0.005),
                          )
        else:
            scales = grzscales
        
    h,w = imgs[0].shape
    rgb = np.zeros((h,w,3), np.float32)
    # Convert to ~ sigmas
    for im,band in zip(imgs, bands):
        plane,scale = scales[band]
        rgb[:,:,plane] = (im / scale).astype(np.float32)
        #print 'rgb: plane', plane, 'range', rgb[:,:,plane].min(), rgb[:,:,plane].max()

    if mnmx is None:
        mn,mx = -3, 10
    else:
        mn,mx = mnmx

    if arcsinh is not None:
        def nlmap(x):
            return np.arcsinh(x * arcsinh) / np.sqrt(arcsinh)
        rgb = nlmap(rgb)
        mn = nlmap(mn)
        mx = nlmap(mx)

    rgb = (rgb - mn) / (mx - mn)

    if desaturate:
        # optionally desaturate pixels that are dominated by a single
        # colour to avoid colourful speckled sky

        RGBim = np.array([rgb[:,:,0],rgb[:,:,1],rgb[:,:,2]])
        a = RGBim.mean(axis=0)
        np.putmask(a, a == 0.0, 1.0)
        acube = np.resize(a,(3,h,w))
        bcube = (RGBim / acube) / 2.5
        mask = np.array(bcube)
        wt = np.max(mask,axis=0)
        np.putmask(wt, wt > 1.0, 1.0)
        wt = 1 - wt
        wt = np.sin(wt*np.pi/2.0)
        temp = RGBim * wt + a*(1-wt) + a*(1-wt)**2 * RGBim
        rgb = np.zeros((h,w,3), np.float32)
        for idx,im in enumerate((temp[0,:,:],temp[1,:,:],temp[2,:,:])):
            rgb[:,:,idx] = im

    clipped = np.clip(rgb, 0., 1.)

    return clipped
 
def run_nsa(nsa_decals,dr='2',random_samp=True,force_fits=False):

    if random_samp:
        N = 101
        galaxies = random.sample(nsa_decals,N)
    else:
        galaxies = nsa_decals

    # Set new parameters

    _scales = dict(g = (2, 0.008), r = (1, 0.014), z = (0, 0.019))
    _mnmx = (-0.5,300)

    good_images = np.ones(len(galaxies),dtype=bool)

    for i,gal in enumerate(galaxies):

        # Check if FITS image already exists
        #volpath = '/Volumes/3TB/gz4/DECaLS'
        volpath = '/Users/willettk/Astronomy/Research/GalaxyZoo/decals'

        fitspath = '%s/fits/nsa' % volpath
        fits_filename = '%s/%s.fits' % (fitspath,gal['IAUNAME'])
        if os.path.exists(fits_filename) == False or force_fits:
            get_skyserver_fits(gal,fitspath,dr,remove_multi_fits=False)

        # Check if new JPEG image already exists
        jpegpath = '%s/jpeg/dr2' % volpath
        jpeg_filename = '%s/%s.jpeg' % (jpegpath,gal['IAUNAME'])
        if os.path.exists(jpeg_filename) == False:
            try:
                img,hdr = fits.getdata(fits_filename,0,header=True)

                badmax = 0.
                for j in range(img.shape[0]):
                    band = img[j,:,:]
                    nbad = (band == 0.).sum() + np.isnan(band).sum()
                    fracbad = nbad / np.prod(band.shape)
                    badmax = max(badmax,fracbad)

                if badmax < 0.2:
                    rgbimg = dstn_rgb((img[0,:,:],img[1,:,:],img[2,:,:]), 'grz', mnmx=_mnmx, arcsinh=1., scales=_scales, desaturate=True)
                    plt.imsave(jpeg_filename, rgbimg, origin='lower')
                else:
                    good_images[i] = False
            except IOError:
                print "Couldn't find image for %s" % fits_filename
                good_images[i] = False

        if not i % 10 and i > 0:
            print '%5i/%5i galaxies completed; %s' % (i,len(galaxies),datetime.datetime.now().strftime("%H:%M:%S"))

    # Write good images to file
    nsa_version = '1_0_0'
    galaxies.write('%s/decals/fits/nsa_v%s_decals_dr%s_goodimgs.fits' % (gzpath,nsa_version,dr),overwrite=True)

    return None

if __name__ == "__main__":

    dr = '2'
    nsa = get_nsa_images(highz=True)
    bricks = get_decals_bricks(dr)
    nsa_decals = run_all_bricks(nsa,bricks,dr,run_to=-1)
    run_nsa(nsa_decals,dr,random_samp=False,force_fits=True)

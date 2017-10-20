import os
import urllib
from multiprocessing import Value, Lock
from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
import progressbar2 as pb
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt

from python.decals_dr2 import dstn_rgb

# analagous to the middle part of decals_dr2.py: given a catalog, download FITS and make JPEG
# Not exactly identical though

widgets = ['Downloads: ', pb.Percentage(), ' ', pb.Bar(marker='0',left='[',right=']'), ' ', pb.ETA()]
cdx = 0
pbar = 0
class Counter(object):
    def __init__(self, initval=0):
        self.val = Value('i', initval)
        self.lock = Lock()
    def change_value(self,v):
        with self.lock:
            self.val.value = v
    def increment(self):
        with self.lock:
            self.val.value += 1
    def value(self):
        with self.lock:
            return self.val.value

min_pixelscale = 0.10
def get_skyserver_fits(gal,fitspath='../fits/nsa',dr='2',remove_multi_fits=False):
    timed_out = False
    good_images = False
    # Download a multi-plane FITS image from the DECaLS skyserver

    # Get FITS
    galname = gal['IAUNAME']
    fits_filename = "{0}/{1}.fits".format(fitspath, galname)
    if os.path.exists(fits_filename) == False:
        params = urllib.urlencode({'ra':gal['RA'],'dec':gal['DEC'],'pixscale':max(min(gal['PETROTH50']*0.04,gal['PETROTH90']*0.02),min_pixelscale),'size':424})
        if dr == '1':
            url = "http://imagine.legacysurvey.org/fits-cutout-decals-dr1?{0}".format(params)
        elif dr == '2':
            url = "http://legacysurvey.org/viewer/fits-cutout-decals-dr2?{0}".format(params)
        try:
            urllib.urlretrieve(url, fits_filename)

            # Write multi-plane FITS images to separate files for each band

            data,hdr = fits.getdata(fits_filename,0,header=True)
            for idx,plane in enumerate('grz'):
                hdr_copy = hdr.copy()
                hdr_copy['NAXIS'] = 2
                hdr_copy['FILTER'] = '{0}       '.format(plane)
                for badfield in ('BANDS','BAND0','BAND1','BAND2','NAXIS3'):
                    hdr_copy.remove(badfield)
                fits.writeto("{0}_{1}.fits".format(fits_filename,plane),data[idx,:,:],hdr_copy,clobber=True)

            if remove_multi_fits:
                os.remove(fits_filename)
            del data,hdr
            timed_out, good_images = makejpeg(gal, fits_filename)
        except:
            timed_out = True
            cdx.increment()
            pbar.update(cdx.value())
    else:
        timed_out, good_images = makejpeg(gal, fits_filename)
    return [timed_out, good_images]

def makejpeg(gal, fits_filename, jpegpath="../jpeg/dr2"):
    good_image = False
    timed_out = False
    _scales = dict(g = (2, 0.008), r = (1, 0.014), z = (0, 0.019))
    _mnmx = (-0.5,300)
    jpeg_filename = '{0}/{1}.jpeg'.format(jpegpath,gal['IAUNAME'])
    if os.path.exists(jpeg_filename) == False:
        if os.path.exists(fits_filename):
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
                    good_image = True
            except:
                timed_out = True
        else:
            timed_out = True
    else:
        # JPEG file already exists, so it's assumed to be a good image
        good_image = True
    cdx.increment()
    pbar.update(cdx.value())
    return [timed_out, good_image]

if __name__ == "__main__":
    dr = '2'
    nsa_version = '1_0_0'
    nsa_decals = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}_after_cuts.fits'.format(nsa_version, dr), 1))
    cdx=Counter(0)
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(nsa_decals))
    pool=ThreadPool(8)
    pbar.start()
    results = pool.map(get_skyserver_fits, nsa_decals)
    # results = map(get_skyserver_fits, nsa_decals)
    pbar.finish()
    pool.close()
    pool.join()
    results = np.array(results)
    timed_out = results[:,0]
    good_images = results[:,1]

    logfile = "../failed_fits_downloads.log"
    flog = open(logfile, 'w')
    print(flog, "\n".join(nsa_decals['IAUNAME'][timed_out]))
    flog.close()

    logfile = "../bad_pix_images.log"
    ilog = open(logfile, 'w')
    print(ilog, "\n".join(nsa_decals['IAUNAME'][~good_images]))
    ilog.close()

    print("\n{0} total galaxies processed".format(len(nsa_decals)))
    print("{0} good images".format(sum(good_images)))
    print("{0} galaxies with bad pixels".format(len(nsa_decals) - sum(good_images)))
    print("{0} galaxies timed out downloading data from Legacy Skyserver".format(sum(timed_out)))

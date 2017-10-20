# Make metadata table for DECaLS images

import os
import warnings

import numpy as np
from astropy import units as u
from astropy.cosmology import WMAP9
from astropy.io import fits
from astropy.table import Table

warnings.simplefilter("ignore",RuntimeWarning)

version = '1_0_0'
nsa_decals = fits.getdata('../fits/nsa_v{0}_decals_dr2_after_cuts_after_failed_removed.fits'.format(version),1)

N = len(nsa_decals)

t = Table()

t['coords.0'] = nsa_decals['RA']
t['coords.1'] = nsa_decals['DEC']

# Calculate absolute size in kpc

size = WMAP9.kpc_proper_per_arcmin(nsa_decals['Z']).to(u.kpc/u.arcsec)*(nsa_decals['PETROTHETA']*u.arcsec)
size[nsa_decals['Z']<0] = -99.*u.kpc

# Calculate absolute and apparent magnitude

absmag_r = nsa_decals['ABSMAG'][:,4].astype(float)
mag = 22.5 - 2.5*np.log10(nsa_decals['NMGY']).astype(float)
mag[~np.isfinite(mag)] = -99.

fluxarr = nsa_decals['PETROFLUX'][:,4].astype(float)

url_stub = "http://www.galaxyzoo.org.s3.amazonaws.com/subjects/decals_dr2"

for imgtype in ('standard','inverted','thumbnail'):
    lst = url_stub + '/{0}/'.format(imgtype) + nsa_decals['IAUNAME'] + '.jpeg'
    t['location.{0}'.format(imgtype)] = lst

t['nsa_id'] = np.char.add('NSA_',nsa_decals['NSAID'].astype(str))

t['metadata.absolute_size'] = size.value
t['metadata.counters.feature'] = np.zeros(N,dtype=int)
t['metadata.counters.smooth'] = np.zeros(N,dtype=int)
t['metadata.counters.star'] = np.zeros(N,dtype=int)

t['metadata.mag.faruv'] = mag[:,0]
t['metadata.mag.nearuv'] = mag[:,1]
t['metadata.mag.u'] = mag[:,2]
t['metadata.mag.g'] = mag[:,3]
t['metadata.mag.r'] = mag[:,4]
t['metadata.mag.i'] = mag[:,5]
t['metadata.mag.z'] = mag[:,6]
t['metadata.mag.abs_r'] = absmag_r

t['metadata.nobs_max.g'] = nsa_decals['nobs_max_g']
t['metadata.nobs_max.r'] = nsa_decals['nobs_max_r']
t['metadata.nobs_max.z'] = nsa_decals['nobs_max_z']

t['metadata.petrorad_50_r'] = nsa_decals['PETROTH50']
t['metadata.petroflux_r'] = fluxarr
t['metadata.petrorad_r'] = nsa_decals['PETROTHETA']
t['metadata.redshift'] = nsa_decals['Z']
t['metadata.retire_at'] = [40] * N
t['survey'] = ['decals'] * N

# Check format and content
#t[:5].show_in_browser()

# Write to table
ftypes = ('csv','fits')
for ft in ftypes:
    fname = '../fits/nsa_decals_dr2_metadata.{0}'.format(ft)
    if os.path.isfile(fname):
        os.remove(fname)
    t.write(fname)

# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import errno
import os
import random
import urllib

import numpy as np
import progressbar as pb
from astropy.io import fits
import astropy.table
from astropy.table import Table
from matplotlib import pyplot as plt

min_pixelscale = 0.10


def get_nsa_catalog(nsa_catalog_loc):
    '''
    Return the loaded NASA-Sloan Atlas catalog

    Args:
        nsa_catalog_loc (str): absolute file path to NSA catalog e.g. [dir]/nsa_v_0_1_2.fits

    Returns:

    '''

    nsa = fits.getdata(nsa_catalog_loc, 1)

    return nsa


def get_decals_bricks(bricks_loc, dr):
    '''

    Return the loaded catalog of DECaLS image bricks
    Include only bricks with images in r, g, and z-bands.
    No constraints on catalog entry for DR2, oddly.

    Args:
        dr (str): load bricks from (dr) data release

    Returns:
        (dict) catalog of RA, Dec edges of DECALS 'brick' tile images in rgz bands
    '''

    bricks_all = fits.getdata(bricks_loc, 1)

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

    elif dr == '3':
        has_g = bricks_all['nexp_g'] > 0
        has_r = bricks_all['nexp_r'] > 0
        has_z = bricks_all['nexp_z'] > 0

        bricks = bricks_all[has_g & has_r & has_z]

    else:
        raise Exception('Data release "{}" not recognised'.format(dr))

    assert len(bricks) > len(bricks_all) * 0.1

    return bricks


def find_matching_brick(gal, bricks):
    '''
    Find the DECaLS brick covering the (RA,dec) coordinates for a galaxy's position
    Galaxy coordinates may be part of multiple bricks

    Args:
        gal (dict): galaxy coordinates in form {RA, DEC}
        bricks (): catalog table of DECALS bricks, listing RA and DE edges

    Returns:
        nmatch (int) total number of bricks containing the galaxy coordinates
        coomatch (int) boolean array of len(bricks), 1 where brick matches galaxy coordinates
    '''

    '''
    Bricks are roughly 0.25 deg x 0.25 deg on each side
    Three sets of coordinates:
    ra,dec = center of brick
    ra1,dec1 = lower right corner of brick
    ra2,dec2 = upper left corner of brick
    '''

    ragal, decgal = gal['RA'], gal['DEC']

    # Find boolean array of bricks that match in RA
    ramatch = (bricks['ra1'] < ragal) & (bricks['ra2'] >= ragal)
    # Find boolean array of bricks that match in dec
    decmatch = (bricks['dec1'] < decgal) & (bricks['dec2'] >= decgal)

    coomatch = (ramatch & decmatch)  # boolean array of bricks that match in both RA and dec
    nmatch = sum(coomatch)

    return nmatch, coomatch


def run_all_bricks(nsa, bricks, dr, nsa_version, run_to=-1):
    '''
    Create a matched catalogue of all NSA sources that have grz imaging in DECaLS

    Galaxies must match the galaxy selection criteria (not yet implemented - 'with_cuts' or 'clean' catalog?)

    Args:
        nsa ():
        bricks ():
        dr (str): DECALS data release version e.g. '2'
        nsa_version (str): NASA Sloan Atlas version e.g. TODO
        run_to (int): nsa galaxies to match. If -1, matches all

    Returns:
        (astropy.Table) of format [{NSA galaxy, NSA details, DECALS brick with that galaxy}]

    '''

    # nsa is dict of arrays of e.g. galaxy ra

    brick_maxdec = max(bricks['dec2'])
    brick_mindec = min(bricks['dec1'])

    # Rough limits of DR1 in RA
    # ralim = ((nsa['RA'] > 15/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 1/24. * 360) | (nsa['RA'] > 21/24. * 360) | (nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 11/24. * 360)
    # declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

    # Make this routine somewhat quicker by first eliminating everything in the NSA catalog
    # outside the observed RA/dec range of the DECaLS bricks.
    # Rough limits of DR2 in RA

    ralim = ((nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 3/24. * 360) | (nsa['RA'] > 21/24. * 360)
    declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

    total_matches = 0
    multi_matches = 0
    decals_indices = np.zeros(len(nsa),
                              dtype=bool)  # binary list of all nsa galaxies, will be 1 if that galaxy is in decals
    bricks_indices = []  # empty list, will contain all coordinates of nsa galaxies matched to (exactly) 1 decals brick

    widgets = ['Matching: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(nsa[:run_to]))
    pbar.start()

    # For every galaxy in the NSA catalog, if it is in DECALS RA/DEC window, find which brick(s) it is in
    for idx, gal in enumerate(nsa[:run_to]):
        if (declim & ralim)[idx]:
            nm, coomatch = find_matching_brick(gal, bricks)
            #  record if one matching brick, or many
            if nm > 0:
                total_matches += 1
                decals_indices[idx] = True  # set nsa filter to True for this galaxy
                bricks_indices.append(coomatch.argmax())  # record the bricks index of the one matching brick
            if nm > 1:
                multi_matches += 1  # count multi-matches but do not consider them as matches
        pbar.update(idx + 1)
    pbar.finish()

    print('{0:6d} total matches between NASA-Sloan Atlas and DECaLS DR{1}'.format(total_matches, dr))
    print('{0:6d} galaxies had matches in more than one brick'.format(multi_matches))

    nsa_table = Table(nsa)  # load NSA catalog into Table (previously accessed as dict-like FitsRecord)
    nsa_decals = nsa_table[decals_indices]  # filter table to only galaxies which are matched to DECALS

    # ought to have that many coordinates in the 'bricks_indices' list
    assert len(nsa_decals) == len(bricks_indices), \
        "Length of nsa_decals ({0}) and bricks_indices ({1}) must match".format(len(nsa_decals), len(bricks_indices))

    matched_bricks = Table(bricks[bricks_indices])

    assert len(nsa_decals) == len(matched_bricks)

    # add the bricks data into the nsa_decals table (manual merge!)
    nsa_decals_bricks = astropy.table.hstack([nsa_decals, matched_bricks])
    assert len(nsa_decals_bricks) == len(matched_bricks)

    # Write to file
    # Check what version of the NSA is being used and set string variable below
    outfilename = '../fits/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, dr)
    nsa_decals.write(outfilename, overwrite=True)

    return nsa_decals


def run_nsa(nsa_decals, dr='2', nsa_version='1_0_0', random_samp=True, force_fits=False):
    '''
    For galaxies with coverage in the DECaLS bricks, download FITS images from cutout service and make JPGs

    Args:
        nsa_decals (astropy.Table): catalog of galaxies in both nsa and decals including RA, Dec, Petro.
        dr (int): data release to download
        nsa_version (?): ?
        random_samp (bool): if True, restrict to a random sample of 101 galaxies only
        force_fits (bool): if True, download FITS even if identically-named file already exists

    Returns:
        None
    '''

    if random_samp:
        N = 101
        galaxies = random.sample(nsa_decals, N)
    else:
        galaxies = nsa_decals

    # Set parameters for RGB image creation
    _scales = dict(
        g=(2, 0.008),
        r=(1, 0.014),
        z=(0, 0.019))

    _mnmx = (-0.5, 300)

    good_images = np.zeros(len(galaxies), dtype=bool)

    # Set paths for the FITS and JPG files

    fitspath = '../fits/nsa'
    make_sure_path_exists(fitspath)

    jpegpath = '../jpeg/dr2'
    make_sure_path_exists(jpegpath)

    # log file, nice trick that I should learn
    logfile = "../failed_fits_downloads.log"
    flog = open(logfile, 'w')
    timed_out = np.zeros(len(galaxies), dtype=bool)

    widgets = ['Downloads: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(galaxies))
    pbar.start()
    for i, gal in enumerate(galaxies):

        # Check if FITS image already exists
        fits_filename = '{0}/{1}.fits'.format(fitspath, gal['IAUNAME'])
        if os.path.exists(fits_filename) is False or force_fits:

            # Download the FITS image
            try:
                get_skyserver_fits(gal, fitspath, dr, remove_multi_fits=False)
            except IOError:
                print("IOError downloading {0}".format(gal['IAUNAME']))
                timed_out[i] = True
                print(flog, gal['IAUNAME'])

        # Check if new JPEG image already exists
        jpeg_filename = '{0}/{1}.jpeg'.format(jpegpath, gal['IAUNAME'])
        # if FITS exists and JPEG does not
        if os.path.exists(jpeg_filename) is False:
            if os.path.exists(fits_filename):
                # generate the JPEG from FITS data using Dustin's RGB image creator
                try:
                    img, hdr = fits.getdata(fits_filename, 0, header=True)

                    badmax = 0.
                    for j in range(img.shape[0]):
                        band = img[j, :, :]
                        nbad = (band == 0.).sum() + np.isnan(band).sum()
                        fracbad = nbad / np.prod(band.shape)
                        badmax = max(badmax, fracbad)

                    if badmax < 0.2:
                        rgbimg = dstn_rgb(
                            (img[0, :, :], img[1, :, :], img[2, :, :]), 'grz',
                            mnmx=_mnmx,
                            arcsinh=1.,
                            scales=_scales,
                            desaturate=True)
                        plt.imsave(jpeg_filename, rgbimg, origin='lower')
                        good_images[i] = True
                except:
                    print("Other download error for {0}".format(gal['IAUNAME']))
                    timed_out[i] = True
                    print(flog, gal['IAUNAME'])
            else:
                print("Could not find FITS file for {0}".format(fits_filename))
        else:
            # JPEG file already exists, so it's assumed to be a good image
            good_images[i] = True
        pbar.update(i + 1)

    pbar.finish()

    # Close logging file for timed-out errors from server
    flog.close()
    galaxies[timed_out].write('../fits/nsa_v{0}_decals_dr{1}_timedout.fits'.format(nsa_version, dr), overwrite=True)

    # Write good images to file
    galaxies[good_images].write('../fits/nsa_v{0}_decals_dr{1}_goodimgs.fits'.format(nsa_version, dr), overwrite=True)

    # Print summary to screen

    print("\n{0} total galaxies processed".format(len(galaxies)))
    print("{0} good images".format(sum(good_images)))
    print("{0} galaxies with bad pixels".format(len(galaxies) - sum(good_images)))
    print("{0} galaxies timed out downloading data from Legacy Skyserver".format(sum(timed_out)))

    return None


def get_skyserver_fits(gal, fitspath, dr='1', remove_multi_fits=True):
    '''
    Download a multi-plane FITS image from the DECaLS skyserver
    Adapt download pixel scale by catalog radius
    Save separate files in each filter and optionally delete multi-band

    Args:
        gal (dict): galaxy coordinates {RA, DEC}
        fitspath (str): location to save files
        dr (int): which data release to download
        remove_multi_fits (bool): if True, remove the multi-band images afterwards

    Returns:
        None
    '''

    # Get FITS
    galname = gal['IAUNAME']
    params = urllib.parse.urlencode({'ra': gal['RA'], 'dec': gal['DEC'],
                               'pixscale': max(min(gal['PETROTH50'] * 0.04, gal['PETROTH90'] * 0.02), min_pixelscale),
                               'size': 424})
    if dr == '1':
        url = "http://imagine.legacysurvey.org/fits-cutout-decals-dr1?{0}".format(params)
    elif dr == '2':
        url = "http://legacysurvey.org/viewer/fits-cutout-decals-dr2?{0}".format(params)

    urllib.request.urlretrieve(url, "{0}/{1}.fits".format(fitspath, galname))

    # Write multi-plane FITS images to separate files for each band
    # Alternatively, could directly specify the desired band?

    # Some kind of checking for corrupted filter images?
    data, hdr = fits.getdata("{0}/{1}.fits".format(fitspath, galname), 0, header=True)
    for idx, plane in enumerate('grz'):
        hdr_copy = hdr.copy()
        hdr_copy['NAXIS'] = 2
        hdr_copy['FILTER'] = '{0}       '.format(plane)
        for badfield in ('BANDS', 'BAND0', 'BAND1', 'BAND2', 'NAXIS3'):
            hdr_copy.remove(badfield)
        fits.writeto("{0}/{1}_{2}.fits".format(fitspath, galname, plane), data[idx, :, :], hdr_copy, overwrite=True)

    if remove_multi_fits:
        os.remove("{0}/{1}.fits".format(fitspath, galname))
    del data, hdr

    return None


def dstn_rgb(imgs, bands, mnmx=None, arcsinh=None, scales=None, desaturate=False):
    '''
    Given a list of image arrays in the given bands, returns a scaled RGB image.

    Args:
        imgs (list): numpy arrays, all the same size, in nanomaggies
        bands (list): strings, eg, ['g','r','z']
        mnmx (min,max), values that will become black/white *after* scaling. Default is (-3,10)):
        arcsinh (bool): if True, use nonlinear scaling (as in SDSS)
        scales ():
        desaturate ():

    Returns:
        (H,W,3) numpy array with values between 0 and 1.

    '''

    bands = ''.join(bands)

    grzscales = dict(
        g=(2, 0.0066),
        r=(1, 0.01385),
        z=(0, 0.025),
    )

    if scales is None:
        if bands == 'grz':
            scales = grzscales
        elif bands == 'urz':
            scales = dict(
                u=(2, 0.0066),
                r=(1, 0.01),
                z=(0, 0.025),
            )
        elif bands == 'gri':
            scales = dict(
                g=(2, 0.002),
                r=(1, 0.004),
                i=(0, 0.005),
            )
        else:
            scales = grzscales

    h, w = imgs[0].shape
    rgb = np.zeros((h, w, 3), np.float32)
    # Convert to ~ sigmas
    for im, band in zip(imgs, bands):
        plane, scale = scales[band]
        rgb[:, :, plane] = (im / scale).astype(np.float32)

    if mnmx is None:
        mn, mx = -3, 10
    else:
        mn, mx = mnmx

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

        RGBim = np.array([rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]])
        a = RGBim.mean(axis=0)
        np.putmask(a, a == 0.0, 1.0)
        acube = np.resize(a, (3, h, w))
        bcube = (RGBim / acube) / 2.5
        mask = np.array(bcube)
        wt = np.max(mask, axis=0)
        np.putmask(wt, wt > 1.0, 1.0)
        wt = 1 - wt
        wt = np.sin(wt*np.pi/2.0)
        temp = RGBim * wt + a*(1-wt) + a*(1-wt)**2 * RGBim
        rgb = np.zeros((h, w, 3), np.float32)
        for idx, im in enumerate((temp[0, :, :], temp[1, :, :], temp[2, :, :])):
            rgb[:, :, idx] = im

    clipped = np.clip(rgb, 0., 1.)

    return clipped


def make_sure_path_exists(path):

    # Check if a local path exists; if not, create it.

    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


if __name__ == "__main__":

    # Run all steps to create the NSA-DECaLS-GZ catalog

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'

    nsa_version = '0_1_2'
    nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

    # TODO should document which files are being used and then rename
    dr = '3'
    if dr == '3':
        # http: // legacysurvey.org / dr3 / files /
        bricks_filename = 'survey-bricks-dr3-with-coordinates.fits'
    elif dr == '2':
        # http: // legacysurvey.org / dr2 / files /
        bricks_filename = 'decals-bricks-dr2.fits'
    bricks_loc = '{}/{}'.format(catalog_dir, bricks_filename)

    nsa = get_nsa_catalog(nsa_catalog_loc)
    bricks = get_decals_bricks(bricks_loc, dr)
    nsa_decals = run_all_bricks(nsa, bricks, dr, nsa_version, run_to=300)

    # nsa_decals = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}_after_cuts.fits'.format(nsa_version, dr), 1))
    # TODO still need to apply cuts - happened externally via Topcat
    # nsa_decals = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, dr), 1))
    # run_nsa(nsa_decals, dr, random_samp=False, force_fits=False)

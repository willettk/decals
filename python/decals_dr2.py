# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import errno
import os
import random

import numpy as np
from astropy.io import fits
import astropy.table
from astropy.table import Table
from tqdm import tqdm

from python.get_images.download_images_threaded import download_images_multithreaded

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

    # For every galaxy in the NSA catalog, if it is in DECALS RA/DEC window, find which brick(s) it is in
    galaxies_enumerated = tqdm(enumerate(nsa[:run_to]), total=len(nsa[:run_to]), unit=' galaxies checked')
    for idx, gal in galaxies_enumerated:
        if (declim & ralim)[idx]:
            nm, coomatch = find_matching_brick(gal, bricks)
            #  record if one matching brick, or many
            if nm > 0:
                total_matches += 1
                decals_indices[idx] = True  # set nsa filter to True for this galaxy
                bricks_indices.append(coomatch.argmax())  # record the bricks index of the one matching brick
            if nm > 1:
                multi_matches += 1  # count multi-matches but do not consider them as matches

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


def run_nsa(nsa_decals, dr='2', nsa_version='1_0_0', random_samp=True, overwrite=False):
    '''
    For galaxies with coverage in the DECaLS bricks, download FITS images from cutout service and make JPGs

    Args:
        nsa_decals (astropy.Table): catalog of galaxies in both nsa and decals including RA, Dec, Petro.
        dr (str): data release to download
        nsa_version (str): version of NSA catalog being used. Defines output catalog filenames.
        random_samp (bool): if True, restrict to a random sample of 101 galaxies only
        overwrite (bool): if True, download FITS and remake JPEG even if identically-named file(s) already exist

    Returns:
        None
    '''

    if random_samp:
        N = 101
        galaxies = random.sample(nsa_decals, N)
    else:
        galaxies = nsa_decals

    fits_dir = '../fits/nsa'
    make_sure_directory_exists(fits_dir)

    jpeg_dir = '../jpeg/dr2'
    make_sure_directory_exists(jpeg_dir)

    timed_out, good_images = download_images_multithreaded(nsa_decals, dr, fits_dir, jpeg_dir, overwrite=overwrite)

    # Write catalogs of time-outs and good images to file
    galaxies[timed_out].write('../fits/nsa_v{0}_decals_dr{1}_timedout.fits'.format(nsa_version, dr), overwrite=True)
    galaxies[good_images].write('../fits/nsa_v{0}_decals_dr{1}_goodimgs.fits'.format(nsa_version, dr), overwrite=True)


def make_sure_directory_exists(path):
    """
    Check if a local directory exists; if not, create it.
    TODO remove? Dynamic paths are not a friend

    Args:
        path (str): directory to check/make

    Returns:
        None
    """

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
    nsa_decals = run_all_bricks(nsa, bricks, dr, nsa_version, run_to=30)

    # nsa_decals = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}_after_cuts.fits'.format(nsa_version, dr), 1))
    # TODO still need to apply cuts - happened externally via Topcat
    # nsa_decals = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, dr), 1))
    run_nsa(nsa_decals, dr, random_samp=False, overwrite=True)

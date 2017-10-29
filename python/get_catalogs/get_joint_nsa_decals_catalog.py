
import numpy as np
from astropy.io import fits
import astropy.table
from astropy.table import Table
from tqdm import tqdm

import matplotlib.pyplot as plt

# TODO should check that galaxy names are unique within each catalog, here and in unit tests


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
        dr (str): load bricks from (data_release) data release

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


def create_joint_catalog(nsa, bricks, dr, nsa_version, run_to=-1):
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

    print(len(nsa), len(bricks))

    brick_maxdec = max(bricks['dec2'])
    brick_mindec = min(bricks['dec1'])

    brick_maxra = max(bricks['ra2'])
    brick_minra = min(bricks['ra2'])

    # Make this routine somewhat quicker by first eliminating everything in the NSA catalog
    # outside the observed RA/dec range of the DECaLS bricks.

    fig, ((ul, ur), (ll, lr)) = plt.subplots(2, 2, sharex=True)
    ul.hist(bricks['dec'])
    ur.hist(nsa['DEC'])
    ll.hist(bricks['ra'])
    lr.hist(nsa['RA'])
    plt.show()


    # Rough limits of DR1
    if dr == '1':
        # ralim = ((nsa['RA'] > 15/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 1/24. * 360) | (nsa['RA'] > 21/24. * 360) | (nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 11/24. * 360)
        # declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

        ralim = (nsa['RA'] >= brick_minra) & (nsa['RA'] <= brick_maxra)
        declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

        # ralim = np.ones(len(nsa), dtype=bool)
        # declim = np.ones(len(nsa), dtype=bool)

    # Rough limits of DR2
    elif dr == '2':
        ralim = ((nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 3/24. * 360) | (nsa['RA'] > 21/24. * 360)
        declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

    # Rough limits of DR2
    elif dr == '3':
        # TODO this is the DR2 values and urgently needs to be updated
        ralim = ((nsa['RA'] > 7/24. * 360) & (nsa['RA'] < 18/24. * 360)) | (nsa['RA'] < 3/24. * 360) | (nsa['RA'] > 21/24. * 360)
        declim = (nsa['DEC'] >= brick_mindec) & (nsa['DEC'] <= brick_maxdec)

    # print(declim.sum(), ralim.sum())

    nsa_in_decals_area = nsa[declim & ralim]
    print(len(nsa_in_decals_area))

    total_matches = 0
    multi_matches = 0
    decals_indices = np.zeros(len(nsa),
                              dtype=bool)  # binary list of all nsa galaxies, will be 1 if that galaxy is in decals
    bricks_indices = []  # empty list, will contain all coordinates of nsa galaxies matched to (exactly) 1 decals brick

    # For every galaxy in the NSA catalog, if it is in DECALS RA/DEC window, find which brick(s) it is in
    galaxies_enumerated = tqdm(
        enumerate(nsa_in_decals_area[:run_to]),
        total=len(nsa_in_decals_area[:run_to]),
        unit=' galaxies checked')
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
        "Length of joint_catalog ({0}) and bricks_indices ({1}) must match".format(len(nsa_decals), len(bricks_indices))

    matched_bricks = Table(bricks[bricks_indices])

    assert len(nsa_decals) == len(matched_bricks)

    # add the bricks data into the joint_catalog table (manual merge!)
    nsa_decals_bricks = astropy.table.hstack([nsa_decals, matched_bricks])
    assert len(nsa_decals_bricks) == len(matched_bricks)

    return nsa_decals


def apply_selection_cuts(input_catalog):
    """
    Select only galaxies with PETROTHETA > 3 and not within 1e-3 of default value (i.e. bad measurement)

    Args:
        catalog (astropy.Table): Galaxy catalog including NSA information

    Returns:
        (astropy.Table) catalog of galaxies matching selection criteria above

    """

    petrotheta_above_3 = input_catalog['PETROTHETA'] > 3

    '''
    NSA catalog’s PETROTHETA calculation sometimes fails to a ‘default’ value that is related to the annulus used to
    measure the Petrosian radius.
    The relation is (found via fitting the observed ‘snap to’ values to the annulus values):
    PETROTHETA_snap_to = 0.997 * PROFTHETA ** 0.998
    Any galaxies with PETROTHETA within 1e-3 of the snap_to value likely has the wrong size.
    '''
    # TODO discuss with Coleman: PROFTHETA is 15 values, which value (or all)?
    #
    # print(input_catalog['PROFTHETA'])
    # proftheta = np.array(input_catalog['PROFTHETA'])
    # snap_values = 0.997 * np.power(input_catalog['PROFTHETA'], 0.998)  # ** syntax behaves unexpectedly with Column
    # snap_tolerance = 1e-3
    # snap_lower_limit = snap_values - snap_tolerance
    # snap_upper_limit = snap_values + snap_tolerance
    #
    # above_snap_lower_limit = input_catalog['PETROTHETA'] > snap_lower_limit
    # below_snap_upper_limit = input_catalog['PETROTHETA'] < snap_upper_limit
    # within_snap_window = above_snap_lower_limit & below_snap_upper_limit
    #
    # selected_catalog = input_catalog[petrotheta_above_3 & ~ within_snap_window]
    # print(len(input_catalog), len(selected_catalog))
    # return selected_catalog

    return input_catalog[petrotheta_above_3]


if __name__ == "__main__":

    # The below will move to main routine
    # Run all steps to create the NSA-DECaLS-GZ catalog

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'

    nsa_version = '0_1_2'
    nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

    dr = '3'
    bricks_filename = 'survey-bricks-dr3-with-coordinates.fits'
    bricks_loc = '{}/{}'.format(catalog_dir, bricks_filename)

    nsa = get_nsa_catalog(nsa_catalog_loc)
    bricks = get_decals_bricks(bricks_loc, dr)

    joint_catalog = create_joint_catalog(nsa, bricks, dr, nsa_version, run_to=-1)
    # Write to file
    joint_catalog_loc = '{0}/nsa_v{1}_decals_dr{2}.fits'.format(catalog_dir, nsa_version, dr)
    joint_catalog.write(joint_catalog_loc, overwrite=True)

    selected_joint_catalog = apply_selection_cuts(joint_catalog)
    # Write to file
    selected_joint_catalog_loc = '{0}/nsa_v{1}_decals_dr{2}_after_cuts.fits'.format(catalog_dir, nsa_version, dr)
    selected_joint_catalog.write(selected_joint_catalog_loc, overwrite=True)

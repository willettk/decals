
from astropy.table import Table

from a_download_decals.get_images.download_images_threaded import get_fits_loc
from new_subjects.find_new_subjects import find_new_catalog_images
from previous_decals_subjects import link_previous_subjects_with_nsa


def get_new_images(joint_catalog, previous_subjects, nsa_catalog):
    """
    Filter joint catalog to galaxies in DR5 that were not previously classified
    Args:
        joint_catalog ():
        previous_subjects ():
        nsa_catalog ():

    Returns:

    """

    # find_new_catalog_images expects two Table catalogs
    # galaxy zoo catalog should include where the FITS would be

    # TODO temporary fix
    catalog_dir = '/Volumes/external/decals/fits/dr5'

    # Is the fits loc column for previous subjects really necessary? Doesn't it either
    # - follow from join with all subjects
    # - not exist?
    previous_subjects['fits_loc'] = [get_fits_loc(catalog_dir, galaxy) for _, galaxy in
                                                previous_subjects.iterrows()]
    # TODO a bit messy - I need to settle on astropy vs. pandas as soon as I get to metadata stage

    #  convert to astropy.Table
    previous_subjects = [row.to_dict() for _, row in previous_subjects.iterrows()]  # neaten
    previous_subjects = Table(previous_subjects)

    # compare DR{} with previous subjects to see what's new
    return find_new_catalog_images(old_catalog=previous_subjects, new_catalog=joint_catalog)

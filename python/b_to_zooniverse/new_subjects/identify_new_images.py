# from astropy.io import fits
from astropy.table import Table

from a_download_decals.get_images.download_images_threaded import get_fits_loc
from new_subjects.find_new_subjects import find_new_catalog_images
from previous_decals_subjects import get_previous_subjects_with_nsa


def get_new_images(joint_catalog, previous_subjects, nsa_catalog):

    # add nsa info to previous gz subjects downloaded from data dump
    previous_galaxy_zoo_with_nsa = get_previous_subjects_with_nsa(previous_subjects, nsa_catalog)

    # find_new_catalog_images expects two Table catalogs
    # galaxy zoo catalog should include where the FITS would be

    # TODO temporary fix
    catalog_dir = '/Volumes/external/decals/fits/dr5'

    # Is the fits loc column for previous subjects really necessary? Doesn't it either
    # - follow from join with all subjects
    # - not exist?
    previous_galaxy_zoo_with_nsa['fits_loc'] = [get_fits_loc(catalog_dir, galaxy) for _, galaxy in
                                                previous_galaxy_zoo_with_nsa.iterrows()]
    # TODO a bit messy - I need to settle on astropy vs. pandas as soon as I get to metadata stage
    previous_galaxy_zoo_with_nsa = [row.to_dict() for _, row in previous_galaxy_zoo_with_nsa.iterrows()]  # neaten
    previous_galaxy_zoo_with_nsa = Table(previous_galaxy_zoo_with_nsa)

    # compare DR{} with previous subjects to see what's new
    return find_new_catalog_images(old_catalog=previous_galaxy_zoo_with_nsa, new_catalog=joint_catalog)
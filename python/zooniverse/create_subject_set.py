
import functools
from multiprocessing.dummy import Pool as ThreadPool
import ast
import pandas as pd
from tqdm import tqdm
import astropy.table

from panoptes_client import Panoptes, Project, SubjectSet, Subject

from make_decals_metadata import get_key_astrophysical_columns

# TODO temporary fix
from get_images.download_images_threaded import get_jpeg_loc


def create_prototype_subject_set(catalog_of_new, calibration_catalog):
    """
    This is the part that changes depending on what I want to upload
    Upload 5000 new galaxies
    Also upload from the calibration set: all bars, all rings, and 500 others
    Args:
        catalog_of_new (?): table of galaxies not previously uploaded as subjects, with nsa and decals info
        calibration_catalog (astropy.Table): table of expertly classified galaxies, with

    Returns:

    """

    # get 1000 new galaxies
    new_galaxies = catalog_of_new[:2500]
    # TODO temporary fix
    jpeg_dir = '/Volumes/external/decals/jpeg/dr5'
    new_galaxies['jpeg_loc'] = [get_jpeg_loc(jpeg_dir, galaxy) for galaxy in new_galaxies]

    # get all calibration bar and ring galaxies
    bar_galaxies = calibration_catalog[calibration_catalog['has_bar']]
    ring_galaxies = calibration_catalog[calibration_catalog['has_ring']]
    other_galaxies = calibration_catalog[~calibration_catalog['has_bar'] & ~calibration_catalog['has_ring']][:500]

    calibration_galaxies = astropy.table.vstack([bar_galaxies, ring_galaxies, other_galaxies])
    print(calibration_galaxies.colnames)

    print('{} bars, {} rings'.format(len(bar_galaxies), len(ring_galaxies)))

    new_galaxies_manifest = create_manifest_from_joint_catalog(new_galaxies)
    calibration_galaxies_manifest = create_manifest_from_calibration_catalog(calibration_galaxies)

    upload_manifest_to_galaxy_zoo('decals_prototype', new_galaxies_manifest)
    upload_manifest_to_galaxy_zoo('decals_calibration', calibration_galaxies_manifest)


def create_manifest_from_calibration_catalog(catalog):
    """
    convert to look like a joint catalog, then upload as normal
    Args:
        catalog ():

    Returns:

    """
    image_columns = ['dr2_jpeg_loc', 'colour_jpeg_loc']
    fake_joint_catalog = convert_calibration_to_joint_catalog(catalog, image_columns)
    print(fake_joint_catalog['jpeg_loc'])

    return create_manifest_from_joint_catalog(fake_joint_catalog)


def convert_calibration_to_joint_catalog(calibration_catalog, image_columns):
    """
    Expand calibration catalog to have one row per image, image under 'jpeg_loc'.
    Record which image processing was used under 'selected_image'.
    Args:
        calibration_catalog (pd.DataFrame): expertly-classified galaxies with multiple calibration images
        image_columns (list): strings of column names with locations of possible images e.g. [dr2_jpeg_loc, etc]

    Returns:
        (pd.DataFrame)
    """
    all_calibration_sets = []
    for image_col in image_columns:
        calibration_set = calibration_catalog.copy()
        calibration_set['jpeg_loc'] = calibration_set[image_col]
        calibration_set['selected_image'] = image_col

        all_calibration_sets.append(calibration_set)

    return astropy.table.vstack(all_calibration_sets)


def create_manifest_from_joint_catalog(catalog):
    """
    Create dict of files and metadata
    Metadata filled with key astro data
    Args:
        catalog ():

    Returns:
        (dict) of form {fileloc: {metadata_col: metadata_value}}
    """

    key_data = get_key_astrophysical_columns(catalog).to_pandas()
    key_data_as_dicts = key_data.apply(lambda x: x.to_dict(), axis=1).values

    jpeg_locs = pd.Series(catalog['jpeg_loc']).values

    manifest = list(zip(jpeg_locs, key_data_as_dicts))

    return manifest


def upload_manifest_to_galaxy_zoo(subject_set_name, manifest, galaxy_zoo_id='5733'):

    # Important - don't commit the password!
    zooniverse_login = read_data_from_txt('zooniverse_login.txt')
    Panoptes.connect(**zooniverse_login)

    galaxy_zoo = Project.find(galaxy_zoo_id)

    subject_set = SubjectSet()

    subject_set.links.project = galaxy_zoo
    subject_set.display_name = subject_set_name

    subject_set.save()

    pbar = tqdm(total=len(manifest), unit=' subjects uploaded')

    save_subject_params = {
        'project': galaxy_zoo,
        'pbar': pbar
    }
    save_subject_partial = functools.partial(save_subject, **save_subject_params)

    pool = ThreadPool(30)
    new_subjects = pool.map(save_subject_partial, manifest)
    pbar.close()
    pool.close()
    pool.join()

    print(new_subjects)
    subject_set.add(new_subjects)


def save_subject(manifest_item, project, pbar=None):
        subject = Subject()

        subject.links.project = project
        # TODO should probably do with named tuple
        subject.add_location(manifest_item[0])
        subject.metadata.update(manifest_item[1])

        subject.save()

        if pbar:
            pbar.update()

        return subject


def read_data_from_txt(file_loc):
    """
    Read and evaluate a python data structure saved as a txt file.
    Args:
        file_loc (str): location of file to read

    Returns:
        data structure contained in file
    """
    with open(file_loc, 'r') as f:
        s = f.read()
        return ast.literal_eval(s)


from make_decals_metadata import get_key_astrophysical_columns

import pandas as pd
from panoptes_client import Panoptes, Project, SubjectSet, Subject
from tqdm import tqdm

import ast
import functools
from multiprocessing.dummy import Pool as ThreadPool

import astropy.table


def create_manifest_from_calibration_catalog(catalog, image_columns):
    """
    convert to look like a joint catalog, then upload as normal
    Args:
        catalog ():

    Returns:

    """
    fake_joint_catalog = convert_calibration_to_joint_catalog(catalog, image_columns)

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

import ast
import functools
from multiprocessing.dummy import Pool as ThreadPool
import warnings

import astropy.table
import pandas as pd
from panoptes_client import Panoptes, Project, SubjectSet, Subject
from tqdm import tqdm

from do_upload.make_decals_metadata import get_key_astrophysical_columns
from settings import zooniverse_login_loc


def create_manifest_from_calibration_catalog(catalog, image_columns):
    """
    Convert manifest to look like a joint catalog (one row per galaxy, png_loc column), then upload as normal
    Args:
        catalog (astropy.Table):
        image_columns (list): strings of column names with locations of possible images e.g. [dr2_png_loc, etc]

    Returns:

    """
    fake_joint_catalog = convert_calibration_to_joint_catalog(catalog, image_columns)

    return create_manifest_from_joint_catalog(fake_joint_catalog)


def convert_calibration_to_joint_catalog(calibration_catalog, image_columns):
    """
    Expand calibration catalog to have one row per image, image under 'png_loc'.
    Record which image processing was used under 'selected_image'.
    Args:
        calibration_catalog (astropy.Table): expertly-classified galaxies with multiple calibration images
        image_columns (list): strings of column names with locations of possible images e.g. [dr2_png_loc, etc]

    Returns:
        (pd.DataFrame) calibration catalog to upload, one galaxy per row, 'png_loc' with selected image
    """
    all_calibration_sets = []
    for image_col in image_columns:
        calibration_set = calibration_catalog.copy()
        calibration_set['png_loc'] = calibration_set[image_col]
        calibration_set['selected_image'] = image_col

        all_calibration_sets.append(calibration_set)

    return astropy.table.vstack(all_calibration_sets)


def create_manifest_from_joint_catalog(catalog):
    """
    Create dict of files and metadata
    Catalog including 'png_loc' and key astro data, one galaxy per row
    Args:
        catalog (astropy.Table): catalog to upload

    Returns:
        (dict) of form {png_loc: img.png, key_data: {metadata_col: metadata_value}}
    """

    key_data = get_key_astrophysical_columns(catalog).to_pandas()
    # calibration catalog can have 'selected image' column
    try:
        key_data['selected_image'] = catalog['selected_image']
    except KeyError:
        pass
    key_data_as_dicts = key_data.apply(lambda x: x.to_dict(), axis=1).values

    png_locs = pd.Series(catalog['png_loc']).values

    data = zip(png_locs, key_data_as_dicts)
    manifest = list(map(lambda x: {'png_loc': x[0], 'key_data': x[1]}, data))

    return manifest


def upload_manifest_to_galaxy_zoo(subject_set_name, manifest, galaxy_zoo_id='5733'):
    """
    Save manifest (set of galaxies with metadata prepared) to Galaxy Zoo

    Args:
        subject_set_name (str): name for subject set
        manifest (list): containing dicts of form {png_loc: img.png, key_data: {metadata_col: metadata_value}}
        galaxy_zoo_id (str): panoptes project id e.g. '5733' for Galaxy Zoo

    Returns:
        None
    """

    if 'TEST' in subject_set_name:
        warnings.warn('Testing mode detected - not uploading!')
        return manifest

    # Important - don't commit the password!
    zooniverse_login = read_data_from_txt(zooniverse_login_loc)
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

    subject_set.add(new_subjects)

    return manifest  # for debugging only


def save_subject(manifest_item, project, pbar=None):
    """

    Add manifest item to project. Note: follow with subject_set.add(subject) to associate with subject set.
    Args:
        manifest_item (dict): of form {png_loc: img.png, key_data: some_data_dict}
        project (str): project to upload subject too e.g. '5773' for Galaxy Zoo
        pbar (tqdm.tqdm): progress bar to update. If None, no bar will display.

    Returns:
        None
    """
    subject = Subject()

    subject.links.project = project
    subject.add_location(manifest_item['png_loc'])
    subject.metadata.update(manifest_item['key_data'])

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

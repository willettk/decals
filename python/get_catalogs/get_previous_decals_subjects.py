import json

import pandas as pd
import astropy.table

from tqdm import tqdm

from python.get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog

def get_galaxy_zoo_decals_catalog(subject_loc):
    subjects = pd.read_csv(subject_loc, nrows=100)
    return subjects


def save_decals_subjects_from_subject_data_dump():
    """
    Extract decals subjects from galaxy zoo subject data dump, save to separate file

    Returns:
        None
    """
    data_dump_loc = '~/Downloads/galaxy_zoo_subjects.csv'
    decals_subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'
    data_df = load_subject_data_dump(data_dump_loc)
    decals_data = data_df[data_df['survey'] == 'decals']
    decals_data.to_csv(decals_subject_loc, index=False)


def load_subject_data_dump(data_dump_loc):
    """
    Convenience function to get galaxy zoo data dump with 'metadata' column expanded

    Args:
        data_dump_loc (str): location of data dump

    Returns:
        (pd.DataFrame) subject data dump with 'metadata' column expanded
    """
    data_dump = pd.read_csv(data_dump_loc, nrows=None)
    data_df = split_json_str_to_columns(data_dump, 'metadata')
    return data_df


def split_json_str_to_columns(input_df, json_column_name):
    """
    Expand Dataframe column of json string into many columns

    Args:
        input_df (pd.DataFrame): dataframe with json str column
        json_column_name (str): json string column name

    Returns:
        (pd.DataFrame) input dataframe with json column expanded into many columns
    """
    json_df = pd.DataFrame(list(input_df[json_column_name].apply(json.loads)))
    del input_df[json_column_name]
    return pd.concat([input_df, json_df], axis=1)


if __name__ == "__main__":

    # The below will move to main routine
    # Run all steps to create the NSA-DECaLS-GZ catalog

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'

    nsa_version = '0_1_2'
    nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

    subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'

    nsa = astropy.table.Table(get_nsa_catalog(nsa_catalog_loc))[['IAUNAME', 'RA', 'DEC']]
    galaxy_zoo = get_galaxy_zoo_decals_catalog(subject_loc)

    # Coordinate catalog has uppercase column names. Rename to lowercase match exposure_catalog.

    for colname in nsa.colnames:
        nsa.rename_column(colname, colname.lower())

    nsa = nsa.to_pandas()

    print(nsa.sample().squeeze())
    print(galaxy_zoo.sample().squeeze())

    # crossmatch on RA and DEC
    # for row_index, galaxy in tqdm(enumerate(galaxy_zoo)):
    # merge on brickname, ra, dec by default

    galaxy_zoo_nsa_catalog = pd.merge(galaxy_zoo, nsa, on=['ra', 'dec'], how='left')

    print(len(galaxy_zoo))
    print(len(nsa))
    print(len(galaxy_zoo_nsa_catalog))
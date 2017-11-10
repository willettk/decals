import json

import pandas as pd
import astropy.table

from tqdm import tqdm

from get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog


def get_galaxy_zoo_decals_catalog(subject_loc):
    subjects = pd.read_csv(subject_loc, nrows=None)
    return subjects


def get_previous_subjects_with_nsa(previous_subjects, nsa):
    # TODO break up each concat step into separate functions
    '''

    Args:
        previous_subjects (pd.DataFrame): decals subjects, extracted from galaxy zoo data dump subjects file
        nsa (str):

    Returns:
        (pd.DataFrame)
    '''

    # galaxy zoo catalog has 'provided_object_id' column with IAU name values. Rename to match exposure catalog.
    previous_subjects = previous_subjects.rename(columns={'provided_image_id': 'iauname',
                                            'dr': 'data_release'})

    nsa_id_columns = ['nsa_id', 'iauname', 'ra', 'dec']
    nsa = nsa[nsa_id_columns].to_pandas()
    nsa.loc[:, 'nsa_id'] = nsa['nsa_id'].astype(int)
    assert all(nsa['iauname'].duplicated() == False)
    assert all(nsa['nsa_id'].duplicated() == False)

    # galaxy_zoo_id_columns = ['nsa_id', 'data_release', 'iauname', 'sdss_id', 'sdss_dr7_id', 'sdss_dr8_id', 'sdss_dr12_objid']
    # print(previous_subjects[galaxy_zoo_id_columns].sample(20))

    # split into DR1 (IAU name as ID) and DR2 (NSA ID as ID)
    galaxy_zoo_dr1 = previous_subjects.dropna(subset=['iauname'])
    galaxy_zoo_dr2 = previous_subjects.dropna(subset=['nsa_id'])
    assert len(galaxy_zoo_dr1) + len(galaxy_zoo_dr2) == len(previous_subjects)

    # add nsa data to dr1
    assert all(galaxy_zoo_dr1['data_release'].isnull())
    galaxy_zoo_dr1.loc[:, 'data_release'] = 'DR1'
    assert all(galaxy_zoo_dr1['nsa_id'].isnull())
    del galaxy_zoo_dr1['nsa_id']  # will be filled by the NSA table in merge below
    assert all(galaxy_zoo_dr1['iauname'].duplicated() == False)
    galaxy_zoo_dr1_with_nsa = pd.merge(galaxy_zoo_dr1, nsa, on='iauname', how='inner')
    assert len(galaxy_zoo_dr1_with_nsa) == len(galaxy_zoo_dr1)

    # add nsa data to dr2
    assert all(galaxy_zoo_dr2['data_release'] == 'DR2')
    assert all(galaxy_zoo_dr2['iauname'].isnull())
    del galaxy_zoo_dr2['iauname']  # will be filled by the NSA table in merge below
    galaxy_zoo_dr2.loc[:, 'nsa_id'] = galaxy_zoo_dr2['nsa_id'].str.lstrip('NSA_').astype(int)
    assert all(galaxy_zoo_dr2['nsa_id'].duplicated() == False)

    galaxy_zoo_dr2_with_nsa = pd.merge(galaxy_zoo_dr2, nsa, on='nsa_id', how='inner')
    assert len(galaxy_zoo_dr2_with_nsa) == len(galaxy_zoo_dr2)

    # restack
    galaxy_zoo_with_nsa = pd.concat([galaxy_zoo_dr1_with_nsa, galaxy_zoo_dr2_with_nsa])

    # note that some galaxies will appear twice, correctly, with dif. zooniverse ids

    return galaxy_zoo_with_nsa.set_index('zooniverse_id', drop=True)


if __name__ == "__main__":

    # The below will move to main routine
    # Run all steps to create the NSA-DECaLS-GZ catalog

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'

    # nsa_version = '0_1_2'
    nsa_version = '1_0_0'  # it's crucial to use 1_0_0 for merging on iauname to work - not clear why!
    nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

    subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'

    nsa = get_nsa_catalog(nsa_catalog_loc)
    galaxy_zoo = get_galaxy_zoo_decals_catalog(subject_loc)

    galaxy_zoo_with_nsa = get_previous_subjects_with_nsa(galaxy_zoo, nsa)

    print(galaxy_zoo_with_nsa.sample(10)[['nsa_id', 'iauname', 'data_release', 'ra', 'dec']])

    galaxy_zoo_with_nsa.to_csv('{}/galaxy_zoo_decals_with_nsa_v{}.csv'.format(catalog_dir, nsa_version))

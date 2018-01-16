
import pandas as pd
import json

from astropy.table import Table


def get_previous_decals_subjects(all_subjects, nsa_v1_0_0):
    """
    Extract decals subjects from galaxy zoo (Ouroborous) subject data dump
    Make the metadata consistent
    Requires loading all subjects, hence slow.

    Args:
        all_subjects (pd.DataFrame): subject data dump from Galaxy Zoo complete dump zip (pinned on Slack)
        nsa_v1_0_0 (astropy.Table): NASA-Sloan Atlas of SDSS galaxies. Must be version 1_0_0 to fix metadata

    Returns:
        (astropy.Table): all decals subjects classified by Galaxy Zoo prior to DR5, as a flat table, joined to NSA
    """
    if any(nsa_v1_0_0['nsa_version'] != '1_0_0'):
        raise Exception('Fatal error: using previous subjects requires NSA catalog version 1_0_0')

    data_df = split_json_str_to_columns(all_subjects, 'metadata')
    decals_df = data_df[data_df['survey'] == 'decals']

    # now we would like to clean up the inconsistent metadata and join the decals subjects with the NSA catalog
    decals_and_nsa = link_previous_subjects_with_nsa(decals_df, nsa_v1_0_0)

    #  convert to astropy.Table
    decals_and_nsa = [row.to_dict() for _, row in decals_and_nsa.iterrows()]
    decals_and_nsa = Table(decals_and_nsa)

    return decals_and_nsa


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


def link_previous_subjects_with_nsa(previous_subjects, nsa):
    '''
    Link previous DECALS subjects with the corresponding entry in the NSA v1_0_0 catalog.
    DR1 and DR2 decals subjects as recorded by Galaxy Zoo have limited and inconsistent information.

    Args:
        previous_subjects (pd.DataFrame): decals subjects, extracted/flattened from galaxy zoo data dump subjects file
        nsa (str): NASA-Sloan Atlas catalog version e.g. '0_1_2'

    Returns:
        (pd.DataFrame): Previous DECALS subjects with consistent information, including the NSA catalog id
    '''

    # galaxy zoo catalog has 'provided_object_id' column with IAU name values. Rename to match exposure catalog.
    previous_subjects = previous_subjects.rename(columns={'provided_image_id': 'iauname',
                                            'dr': 'data_release'})

    nsa_id_columns = ['nsa_id', 'iauname', 'ra', 'dec']
    nsa = nsa[nsa_id_columns].to_pandas()
    nsa.loc[:, 'nsa_id'] = nsa['nsa_id'].astype(int)
    assert all(nsa['iauname'].duplicated() == False)
    assert all(nsa['nsa_id'].duplicated() == False)

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

    # note that some galaxies will appear twice, correctly, with dif. b_to_zooniverse ids
    return galaxy_zoo_with_nsa.set_index('zooniverse_id', drop=True)

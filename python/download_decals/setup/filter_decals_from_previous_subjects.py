import json

import pandas as pd


def save_decals_subjects_from_subject_data_dump(all_subjects_loc, decals_subject_loc):
    """
    Extract decals subjects from galaxy zoo subject data dump, save to separate file

    Returns:
        None
    """
    data_dump = pd.read_csv(all_subjects_loc)
    data_df = split_json_str_to_columns(data_dump, 'metadata')
    decals_data = data_df[data_df['survey'] == 'decals']
    decals_data.to_csv(decals_subject_loc, index=False)


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

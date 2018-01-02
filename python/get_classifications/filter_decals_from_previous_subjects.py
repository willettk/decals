import json

import pandas as pd


def get_decals_subjects_from_all_subjects(all_subjects):
    """
    Extract decals subjects from galaxy zoo subject data dump, save to separate file

    Returns:
        None
    """

    data_df = split_json_str_to_columns(all_subjects, 'metadata')
    return data_df[data_df['survey'] == 'decals']



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

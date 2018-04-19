import json

import pandas as pd

from c_from_zooniverse.decals_schema import TASK_ANSWERS


def classifications_to_table(classifications):
    """
    Panoptes classifications are exported as [user, subject, {response to T0, response to T1, etc}]
    Flatten this into [user, subject, response]
    Args:
        classifications (pd.DataFrame): Panoptes classifications export: [user, subject, {response to T0, etc}]

    Returns:
        (pd.DataFrame): Classifications flattened into simple table of [user, subject, response_to_any_question]
    """
    all_flat_data = []
    # each classification is by one user for one subject, with a json of all responses
    for row_n, classification in classifications.iterrows():
        annotations = json.loads(classification['annotations'])  # get the user responses [(task, value)] for a subject
        clean_annotations = remove_markdown(annotations)  # remove annoying markdown from each response
        flat_df = pd.DataFrame(clean_annotations)  # put each [task, value] pair on a dataframe row
        flat_df['user_id'] = classification['user_id']  # record the (same) user and subject on every row
        flat_df['subject_id'] = classification['subject_ids']
        # copy more columns as needed e.g. time, ip
        all_flat_data.append(flat_df)
    # stick together the responses of all users and all subjects
    return pd.concat(all_flat_data).reset_index(drop=True)  # concat messes with the index


def remove_markdown(annotation):
    """
    Each response value includes the markdown text used to display the image. Remove this.
    e.g. Convert '![img.png](https://panoptes-uploads.etc.png) Features or Disk' to 'Features or Disk'
    Args:
        annotation (list): user response pairs of
            'task' e.g. 'T0' and 'value' e.g. ![img.png](url.png) 'Features or Disk'

    Returns:
        (list): user response pairs without markdown
    """
    for task_response in annotation:
        value_str = task_response['value']
        value_str = value_str.split(')')[-1]
        value_str = value_str.lstrip(' ').rstrip(' ').lower()
        task_response['value'] = value_str
    return annotation


def aggregate_classifications(df):
    """
    Reduce table of all responses by all users for one subject to
    a dict of form {task: {answer: [each user_id with that answer]}
    Args:
        df (pd.DataFrame): all responses by all users for one subject

    Returns:
        (dict): of form {task: {answer: [each user_id with that answer]}
    """
    data = {}
    for task in df['task'].unique():
        data[task] = {}
        for answer in TASK_ANSWERS[task]:
            matching_answers = df[(df['task'] == task) & (df['value'] == answer)]

            data[task].update(
                {answer: list(matching_answers['user_id'])}
            )

    # TODO could have each column as e.g. T0_Smooth, rather than nested dict?
    return data


if __name__ == '__main__':

    panoptes_dir = '/data/galaxy_zoo/decals/panoptes'  # location of Panoptes data exports
    classifications_loc = '{}/{}'.format(panoptes_dir, 'decals-workflow-classifications.csv')  # classification export
    subjects_loc = '{}/{}'.format(panoptes_dir, 'galaxy-zoo-subjects.csv')  # subjects export

    nested_classifications = pd.read_csv(classifications_loc)
    subjects = pd.read_csv(subjects_loc)

    classifications_table = classifications_to_table(nested_classifications)
    aggregated_classifications = classifications_table.groupby('subject_id').apply(aggregate_classifications)

    aggregated_classifications.to_csv('temp.csv')

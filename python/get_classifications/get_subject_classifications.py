
import json
import pandas as pd
import numpy as np

TASK_ANSWERS = {
    'T0':  ['smooth', 'features or disk', 'star or artefact'],
    'T1': ['completely round', 'inbetween', 'cigar-shaped'],
    'T2': ['yes', 'no'],
    'T3': ['rounded', 'boxy', 'no bulge'],
    'T4': ['bar', 'no bar'],
    'T5': ['spiral', 'no spiral'],
    'T6': ['tight', 'medium', 'loose'],
    'T7': ['1', '2', '3', '4', 'more than 4', "can't tell"],
    'T8': ['no bulge', 'just noticeable', 'obvious', 'dominant'],
    'T9': ['yes', 'no'],
    'T10': ['ring', 'lens or arc', 'disturbed', 'irregular', 'other', 'merger', 'dust lane']
}


class Subject():

    def __init__(self):
        pass


class Classification():

    def __init__(self):
        pass


def get_task_name_from_id(task_id):
    """

    Args:
        task_id ():

    Returns:

    """
    mapping = {
        'T0': 'features',
        'T1': 'edge',
        'T2': 'has_bar',
        'T3': 'has_spiral',
        'T4': 'prominent_bulge',
        'T5': 'merging',
        'T6': 'odd',
        'T7': 'rounded',
        'T8': 'has_bulge',
        'T9': 'tight_arms',
        'T10': 'count_arms'
    }
    return task_id[mapping]


def create_subject_objects(subject_data):
    pass


def unpack_json_to_columns(df, json_column_name):
    pass


def get_clean_annotations(annotation):
    """
    Each response value includes the markdown text used to display the image. Remove this.
    Args:
        annotation ():

    Returns:

    """
    for task_response in annotation:
        value_str = task_response['value']
        value_str = value_str.split(')')[-1]
        value_str = value_str.lstrip(' ').rstrip(' ').lower()
        task_response['value'] = value_str
    return annotation


def create_response(classification):
    annotations = json.loads(classification['annotations'])
    clean_annotations = get_clean_annotations(annotations)
    return {
        'annotations': clean_annotations,
        'user_id': classification['user_id']
    }


def aggregate_responses(responses):
    """

    Args:
        responses (list): in form { user_id: 'user_a', annotations: [{task: 'T0', value: 'smooth'}] }

    Returns:
        dict of form {T0: {responses: {'smooth': ['user_a', 'user_b']}

    """

    aggregated_responses = {}
    tasks = set([response['task'] for response in responses])
    for task in tasks:

        # get all the responses to that task
        task_responses = []
        for response in responses:
            for annotation in response['annotations']:
                if annotation['task'] == task:
                    task_responses.append({
                        'user_id': response['user_id'],
                        'response': annotation['value']
                    })

        # get all the possible answers to that task
        answer_types = set([task_response['response'] for task_response in task_responses])

        # list users who gave each answer

        # aggregated_responses[task] = {'responses': }


def flatten_classifications(classifications):
    all_flat_data = []
    for row_n, classification in classifications.iterrows():
        annotations = json.loads(classification['annotations'])
        clean_annotations = get_clean_annotations(annotations)
        flat_df = pd.DataFrame(clean_annotations)
        flat_df['user_id'] = classification['user_id']
        flat_df['subject_id'] = classification['subject_ids']
        # copy more columns as needed
        all_flat_data.append(flat_df)
    return pd.concat(all_flat_data)


def get_aggregate_classifications(flat_classifications):
    aggregated = flat_classifications.groupby('subject_id').apply(aggregate_classifications)

    return aggregated


def aggregate_classifications(df):
    data = {}
    for task in df['task'].unique():
        data[task] = {}
        for answer in TASK_ANSWERS[task]:
            matching_answers = df[(df['task'] == task) & (df['value'] == answer)]

            data[task].update(
                {answer: list(matching_answers['user_id'])}
            )

    # data = pd.DataFrame(data)
    return data


if __name__ == '__main__':

    panoptes_dir = '/data/galaxy_zoo/decals/panoptes'
    classifications_loc = '{}/{}'.format(panoptes_dir, 'decals-workflow-classifications.csv')
    subjects_loc = '{}/{}'.format(panoptes_dir, 'galaxy-zoo-subjects.csv')

    classifications = pd.read_csv(classifications_loc)
    subjects = pd.read_csv(subjects_loc)

    flat = flatten_classifications(classifications)
    print(flat)
    aggregated = get_aggregate_classifications(flat)
    aggregated.to_csv('temp.csv')

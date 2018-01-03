import pytest

from get_subject_classifications import *


"""
The 'annotations' of one user on one subject are returned by Panoptes like the dicts below, stored as a JSON
These need to be flattened into a huge table of [user, subject, task, value] for every user, subject, task and value
"""


@pytest.fixture()
def raw_annotations_a():
    return [
        {"task": "T0",  # only in annotation a
         "task_label": "Is the galaxy simply smooth and rounded, with no sign of a disk?",
         "value": "![feature_or_disk.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/f353f2f1-a47e-439d-b9ca-020199162a79.png) Features or Disk"},

        {"task": "T2",  # in both a and b, users agree
         "task_label": "Could this be a disk viewed edge-on?",
         "value":"![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No"},

        {"task": "T9",  # in both a and b, users disagree
         "task_label": "Is there anything odd?",
         "value": "![yes.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/503a6354-7f72-4899-b620-4399dbd5cf93.png) Yes"},
    ]


@pytest.fixture()
def raw_annotations_b():
    return [
        {"task": "T1",  # only in annotation b
         "task_label": "How rounded is it?",
         "value": "![rounded.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/94412557-f564-40b9-9423-1d2e47cb1104.png) Completely round"},

        {"task": "T2",  # in both a and b, users agree
         "task_label": "Could this be a disk viewed edge-on?",
         "value": "![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No"},

        {"task": "T9",  # in both a and b, users disagree
         "task_label": "Is there anything odd?",
         "value": "![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No"}
    ]


# pack each annotation up as JSON to appear just like the Panoptes classifications export
@pytest.fixture()
def raw_classifications(raw_annotations_a, raw_annotations_b):
    # TODO multiple subject ids
    return pd.DataFrame([
        {'user_id': 'a',
         'annotations': json.dumps(raw_annotations_a),
         'subject_ids': 's1'},
        {'user_id': 'b',
         'annotations': json.dumps(raw_annotations_b),
         'subject_ids': 's1'}
    ])


# the annotations should be flattened into a table like this
# of form [user, subject, task, value] for every user, subject, task and value
@pytest.fixture()
def flat_classifications():
    # TODO multiple subject ids
    return pd.DataFrame([

        # expected responses from user a
        {'task': 'T0',
         'value': 'features or disk',
         'user_id': 'a',
         'subject_id': 's1'},
        {'task': 'T2',
         'value': 'no',
         'user_id': 'a',
         'subject_id': 's1'},
        {'task': 'T9',
         'value': 'yes',
         'user_id': 'a',
         'subject_id': 's1'},

        # expected responses from user b
        {'task': 'T1',
         'value': 'completely round',
         'user_id': 'b',
         'subject_id': 's1'},
        {'task': 'T2',
         'value': 'no',
         'user_id': 'b',
         'subject_id': 's1'},
        {'task': 'T9',
         'value': 'no',
         'user_id': 'b',
         'subject_id': 's1'}
    ])


#  make sure extraneous markdown is correctly removed from the 'value' field of each user response
def test_remove_markdown(raw_annotations_a):
    cleaned = remove_markdown(raw_annotations_a)
    assert cleaned[0]['value'] == 'features or disk'
    assert cleaned[1]['value'] == 'no'
    assert cleaned[2]['value'] == 'yes'


def test_flatten_classifications(raw_classifications, flat_classifications):
    classifications = classifications_to_table(raw_classifications)
    print(classifications[['subject_id', 'task', 'user_id', 'value']])
    print('\n')
    print(flat_classifications[['subject_id', 'task', 'user_id', 'value']])
    print('\n')
    assert classifications[['subject_id', 'task', 'user_id', 'value']].equals(
        flat_classifications[['subject_id', 'task', 'user_id', 'value']])


"""
The flat table needs to be aggregated into a single (perhaps nested) value for each subject
This will eventually be used for e.g. finding spirals, making ML labels, etc.
Not much point testing this carefully now - I'm not sure how it should be aggregated!
Code below is for the previous version, where get_aggregate_classifications is the groupby
"""


# def test_get_aggregated_classifications(flat_classifications):
#     aggregated = get_aggregate_classifications(flat_classifications)
#     print(aggregated.iloc[0])
#
#
# def test_aggregate_responses(responses):
#     aggregated = aggregate_responses(responses)
#     expected = [
#         {'task': 'T0',
#          'responses': [
#              {'response': 'Features or Disk',
#               'user_id': ['a', 'b']}
#          ]
#          },
#
#         {'task': 'T2',
#          'responses': [
#              {'response': 'No',
#               'user_id': ['a']}
#          ]
#          },
#
#         {'task': 'T4',
#          'responses': [
#              {'response': 'No bar',
#               'user_id': ['a']},
#              {'response': 'Bar',
#               'user_id': ['b']}
#          ]
#          }
#     ]
#     assert aggregated == expected
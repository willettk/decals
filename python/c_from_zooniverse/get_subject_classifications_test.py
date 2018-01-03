import pytest

from get_subject_classifications import *

#
# @pytest.fixture()
# def response_a():
#     return {'user_id': 'a',
#             'annotations': [
#                                {'task': 'T0',
#                                 'value': 'features or disk'},
#
#                                {'task': 'T2',
#                                 'value': 'no'},
#
#                                {'task': 'T4',
#                                 'value': 'no bar'}
#                            ]
#             }
#
#
# @pytest.fixture()
# def response_b():
#     return {'user_id': 'b',
#             'annotations': [
#                                {'task': 'T0',
#                                 'value': 'features or disk'},
#
#                                {'task': 'T4',
#                                 'value': 'bar'}
#                            ]
#             }


@pytest.fixture()
def responses(classification_a, classification_b):
    return [classification_a, classification_b]


@pytest.fixture()
def raw_annotations_a():
    return [
        {"task":"T0",
         "task_label": "Is the galaxy simply smooth and rounded, with no sign of a disk?",
         "value": "![feature_or_disk.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/f353f2f1-a47e-439d-b9ca-020199162a79.png) Features or Disk"},
        {"task":"T2",
         "task_label":"Could this be a disk viewed edge-on?",
         "value":"![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No"},
        {"task":"T4",
         "task_label":"Is there a sign of a bar feature through the centre of the galaxy?",
         "value":"![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No bar"},
        {"task":"T5",
         "task_label":"Is there any sign of a spiral arm pattern?",
         "value":"![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No spiral"},
        {"task":"T8",
         "task_label":"How prominent is the central bulge, compared with the rest of the galaxy?",
         "value":"![bulge_dominate.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/fcafa04e-f9ad-46f2-8fa2-e92591e3ac5c.png) Dominant"},
        {"task":"T9",
         "task_label":"Is there anything odd?",
         "value":"![yes.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/503a6354-7f72-4899-b620-4399dbd5cf93.png) Yes"},
        {"task":"T10",
         "task_label":"Is the odd feature a ring, or is the galaxy disturbed or irregular?",
         "value":"Dust lane"}
    ]


@pytest.fixture()
def raw_annotations_b():
    return [
        {"task": "T0",
         "task_label": "Is the galaxy simply smooth and rounded, with no sign of a disk?",
          "value": "![rounded.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/94412557-f564-40b9-9423-1d2e47cb1104.png) Smooth"},
         {"task": "T1",
          "task_label": "How rounded is it?",
          "value": "![rounded.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/94412557-f564-40b9-9423-1d2e47cb1104.png) Completely round"},
         {"task": "T9",
          "task_label": "Is there anything odd?",
          "value": "![no.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/096879e1-12ae-4df8-abb8-d4a93bc7797f.png) No"}
    ]


@pytest.fixture()
def flat_classifications():
    # TODO multiple subject ids
    return pd.DataFrame([
        {'task': 'T0',
         'value': 'features or disk',
         'user_id': 'a',
         'subject_id': 's1'},
        {'task': 'T2',
         'value': 'no',
         'user_id': 'a',
         'subject_id': 's1'},
        {'task': 'T4',
         'value': 'no bar',
         'user_id': 'a',
         'subject_id': 's1'},
        {'task': 'T0',
         'value': 'features or disk',
         'user_id': 'b',
         'subject_id': 's1'},
        {'task': 'T4',
         'value': 'bar',
         'user_id': 'b',
         'subject_id': 's1'}
    ])
#
# @pytest.fixture()
# def raw_classification_a(raw_annotations_a):
#     return pd.Series({
#         'user_id': 'a',
#         'annotations': json.dumps(raw_annotations_a)
#     })
#
#
# @pytest.fixture()
# def raw_classification_b(raw_annotations_b):
#     return pd.Series({
#         'user_id': 'b',
#         'annotations': json.dumps(raw_annotations_b)
#     })


@pytest.fixture()
def raw_classifications(raw_annotations_a, raw_annotations_b):
    return pd.DataFrame([
        {'user_id': 'a',
         'annotations': json.dumps(raw_annotations_a),
         'subject_ids': 's1'},
        {'user_id': 'b',
         'annotations': json.dumps(raw_annotations_b),
         'subject_ids': 's1'}
    ])


def test_clean_annotation(raw_annotations_a):
    cleaned = remove_markdown(raw_annotations_a)
    assert cleaned[0]['value'] == 'features or disk'
    assert cleaned[1]['value'] == 'no'
    assert cleaned[2]['value'] == 'no bar'
    assert cleaned[3]['value'] == 'no spiral'
    assert cleaned[4]['value'] == 'dominant'
    assert cleaned[5]['value'] == 'yes'
    assert cleaned[6]['value'] == 'dust lane'

#
# def test_create_response(raw_classification_a):
#     response = create_response(raw_classification_a)
#     assert response['annotations'][0]['value'] == 'features or disk'
#     assert response['annotations'][6]['value'] == 'dust lane'
#     assert response['user_id'] == 'a'


def test_flatten_classifications(raw_classifications):
    flat_classifications = classifications_to_table(raw_classifications)
#     TODO assert statements

#
# def test_get_aggregated_classifications(flat_classifications):
#     aggregated = get_aggregate_classifications(flat_classifications)
#     print(aggregated.iloc[0])

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
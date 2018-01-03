#### What's the meaning of each Galaxy Zoo data export?


##### Classifications

{workflow}-workflow-classifications.csv lists all the responses made on a given workflow.
{project}-classifications does the same for all workflows within the project.

Each row is the decision tree responses of one user to one subject:
- classification id
- user_name (e.g. klmasters), user_id, user_ip
- workflow_id, workflow_name (e.g. decals_workflow), workflow_version
- created_at
- metadata. Details on the user session e.g. browser, screen size, timezone. NOT manifest data.
- annotations. This is the classification answers as JSON.
- subject_data: the manifest data on that subject
- subject_ids: is the non-plural subject_id?


##### Workflows

{project}-workflows.csv describes the workflow questions.
Each row is the questions for one workflow version.


##### Subjects

{project}-subjects.csv lists the subjects in a project. It doesn't aggregate classifications.
Each row is a subject in the project:
- subject and subject set id
- total classification counts
- status, retirement reason
- manifest metadata (e.g. RA, DEC, etc) as JSON

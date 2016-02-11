from pymongo import MongoClient

client = MongoClient('localhost', 27017)

gz = client['galaxy_zoo'] 

subjects = gz['galaxy_zoo_subjects']
classifications = gz['galaxy_zoo_classifications']
users = gz['galaxy_zoo_users']

ouroboros = client['ouroboros']
discussions = ouroboros['discussions']

tags = ['wrongsize','wrong_size','imagesizewrong']
tags = ['decals_red_artifact',]
wstags = discussions.find({'comments.tags':{'$in':tags},'focus.type':'GalaxyZooSubject'})

ids = []
for tag in wstags:
    ids.append(tag['focus']['_id'])

surveys = []

for id in ids:
    s = subjects.find_one({'zooniverse_id':id})
    survey = s['metadata']['survey']
    surveys.append(survey)

    if survey == 'decals':
        print 'http://talk.galaxyzoo.org/#/subjects/%s' % id

from collections import Counter

print Counter(surveys)



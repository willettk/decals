from pymongo import MongoClient

client = MongoClient()
db = client.galaxy_zoo
subjects = db.subjects

# cursor = subjects.find({'classification_count': 40})
cursor = subjects.find({'metadata.survey': 'decals'})
# cursor = subjects.find({'metadata.dr': 'DR2'})

for document in cursor[:10]:
    print(document)
    # metadata = document['metadata']
    # print(metadata)
    # print(metadata['provided_image_id'])
    # print(metadata['dr'])
    # print(document['location']['standard'])

import pymongo

myclient = pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["mp"]

mycol = mydb["mp_last_task_id_for_new"]
mydict = { "last_idx": 1000000 } # use 2000000 for korea server
x = mycol.insert_one(mydict)

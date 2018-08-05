import pickle

with open("new_CMS_likes_partial.pkl", 'rb') as pkl_file: 
    CMSlikes = pickle.load(pkl_file)

print(CMSlikes)
#for i in range(1000):
#    print(CMSlikes[i])

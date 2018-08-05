import pickle

with open("new_CMS_likes.pkl", 'rb') as pkl_file: 
    CMSlikes = pickle.load(pkl_file)

print(CMSlikes[CMSlikes!=0])
#for i in range(1000):
#    print(CMSlikes[i])

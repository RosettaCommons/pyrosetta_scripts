from collections import defaultdict
import time,os
import pickle
#a={1:set([(1,2),(3,4)])}
#b={1:set([(5,6),(7,9)])}
t0=time.time()
hash_name_list=[]
for hash_name in os.listdir(os.getcwd()):
    if hash_name.startswith('PP2_hashtable'):
        hash_name_list.append(hash_name)
print hash_name_list

hashtable=pickle.load(open(hash_name_list[0],"rb"))
for i in range(1,len(hash_name_list)):
    print hash_name_list[i]
    hashtable_2=pickle.load(open(hash_name_list[i],"rb"))
    for k,v in hashtable_2.items():
        hashtable[k]=hashtable[k] | v
    print time.time()-t0

pickle.dump(hashtable,open('PP2_hashtable_all_5','wb'))
#hashtable_2=pickle.load(open('PP2_hashtable_3_2',"rb"))
#print time.time()-t0
#c=[a,b]
#print hashtable_2
#c=[hashtable,hashtable_2]
#dic = defaultdict(set)
#for k,v in hashtable_2.items():
#   hashtable[k]=hashtable[k] | v
    #print k
    #print v

#print dic
#pickle.dump(hashtable,open('PP2_hashtable_all_3','wb'))


#!/home/sheffler/venv/david/bin/python
#written by yihang
import pyrosetta
from math import *
from xyzMath import *

def calc_distance(v1,v2):
    dis = sqrt((v1-v2)*(v1-v2))
    return dis

def get_surface_unit(surface):
    list = []
    for i in range(1,surface.size()+1):
        print surface.residue(i).annotated_name()
        if surface.residue(i).annotated_name()[-2] != 'L':
            continue
        num = 0
        vec_CA = Vec(surface.residue(i).xyz('Ca2p'))
    #print vec_CA
        for j in range(1,surface.size()+1):
            if surface.residue(j).annotated_name()[-2] == 'L':
                continue
            vec_CO = Vec(surface.residue(j).xyz('C1'))
        #print vec_CO
            if i == 289 and j == 295 :
                print calc_distance(vec_CA,vec_CO)
            if calc_distance(vec_CA,vec_CO) >2.8 and calc_distance(vec_CA,vec_CO) <3.0:
                num = num + 1
        if num == 3:
            list.append(int(i))
        

    print list
    print len(list)
    unit_list = []
    for i in range(len(list)):
        CA1_vec = Vec(surface.residue(list[i]).xyz('Ca2p'))
        for j in range(i+1,len(list)):
            CA2_vec = Vec(surface.residue(list[j]).xyz('Ca2p'))
            if calc_distance(CA1_vec,CA2_vec) > 4.9 and calc_distance(CA1_vec,CA2_vec) < 5.1 :
                #print i
                #print j
                vac_temp = []
                num = 0
                for k in range(1,surface.size()+1):
                    if surface.residue(k).annotated_name()[-2] == 'L':
                        continue
                    vec_CO = Vec(surface.residue(k).xyz('C1'))
                    if (calc_distance(CA1_vec,vec_CO) >2.8 and calc_distance(CA1_vec,vec_CO) <3.0 and calc_distance(CA2_vec,vec_CO) >5.7 and calc_distance(CA2_vec,vec_CO) <5.9) or (calc_distance(CA2_vec,vec_CO) >2.8 and calc_distance(CA2_vec,vec_CO) <3.0 and calc_distance(CA1_vec,vec_CO) >5.7 and calc_distance(CA1_vec,vec_CO) <5.9):
                        #print 'yesyes'
                        num = num +1
                        vac_temp.append(k)
                if num ==2:
                    for m in range(0,2):
                        vec_CO = Vec(surface.residue(vac_temp[m]).xyz('C1'))
                        if calc_distance(CA1_vec,vec_CO) >2.8 and calc_distance(CA1_vec,vec_CO) <3.0:
                            unit_tuple = [list[i],list[j],vac_temp[m]]
                            vac_temp.remove(vac_temp[m])
                            unit_tuple.append(vac_temp[0])
                            break
                    unit_list.append(tuple(unit_tuple))
    print unit_list
    print len(unit_list)
    return unit_list


'''
pyrosetta.init("-extra_res_fa CAL.params CO3.params")
surface = pyrosetta.rosetta.core.import_pose.pose_from_file("A011_6by6by8_CleaveSurface_CAL_CO3.pdb")
bb = get_surface_unit(surface)
print bb
'''

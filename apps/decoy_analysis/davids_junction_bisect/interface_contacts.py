from pyrosetta import rosetta,init
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.scoring.dssp import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.core.select.residue_selector import *
import glob
import os


def count_contacts_accross_junction(pose,resN):
    ss=Dssp(pose).get_dssp_secstruct()
    if ss[resN] != 'H':
        print('Warning: junction residue not helix:  %s'%resN)
        return -1
    in_helix,before_helix,after_helix=identify_helical_segments(ss,resN)
    before_contact_res=get_contacts(in_helix,before_helix, after_helix, pose)
    after_contact_res=get_contacts(in_helix,after_helix,before_helix,pose)
    #print(resN,before_contact_res,after_contact_res)
    return len(before_contact_res)+len(after_contact_res)

def get_contacts(helix,set1,set2,pose):
#this is super wastefull as contacts are computed with whole pose.  need to make subpose over just residues in set, or limit neighborhood selector to this region. but small cost in overall scheme of things
    res_selector=ResidueIndexSelector()
    for index in helix:
        res_selector.append_index(index)
    for index in set1:
        res_selector.append_index(index)

    res_indices=res_selector.apply(pose)
    nb_selector= NeighborhoodResidueSelector(res_indices, 8, False)
    nb_indices=nb_selector.apply(pose)
    contact_res=  [ index for index in range(1, len(nb_indices) + 1) if nb_indices[index] ]
    nearby_contact_res=set(contact_res).intersection(set(set2))
    return nearby_contact_res

def identify_helical_segments(ss,resN):
#identify residues in same helix
    in_helix=[]
    resi=resN
    resT=len(ss)
    while ss[resi] == 'H' and resi > 0:
        in_helix.append(resi+1)
        resi=resi-1

    H_begin=resi

    resi=resN
    while ss[resi] == 'H' and resi < resT:
        in_helix.append(resi+1)
        resi=resi+1
    H_end=resi

# identify residues in preceding two helices
    before_helix=[]
    h_index=0
    in_H=False
    for i in range(H_begin-1,0,-1):
        if ss[i] == 'H':
            before_helix.append(i)
            if not in_H:
                h_index=h_index+1
                in_H=True
                if h_index == 3: break
        else:
            in_H = False

# identify residues in following two helices
    after_helix=[]
    h_index=0
    in_H=False
    for i in range(H_end+1,resT):
        if ss[i] == 'H':
            after_helix.append(i)
            if not in_H:
                h_index=h_index+1
                in_H=True
                if h_index == 3: break
        else:
            in_H = False
    return in_helix,before_helix, after_helix

if __name__ == '__main__':
   #init()
   #ss=['L' for i in range(100)]
   #starts=[10,30,50,70,90]
   #ends=[20,40,60,80,100]
   #for i in range(len(starts)):
   #    for j in range(starts[i],ends[i]):
   #        ss[j]='H'
   #print(ss)
   #in_helix,before_helix,after_helix=identify_helical_segments(ss,52)
   #print(in_helix)
   #print(before_helix)
   #print(after_helix)
   #p=rosetta.core.import_pose.pose_from_file('test.pdb')
   #for ires in range(1,180,5):
   #    nc=count_contacts_accross_junction(p,ires)
   #    print('res %s contacts %s'%(ires,nc))
    pdbs = glob.glob("junction_pdbs/*.pdb")
    output = open("tmp_data.txt","w")
    for pdb in pdbs:
        p=rosetta.core.import_pose.pose_from_file(pdb)
        pdb_name=os.path.basename(pdb)[0:-4]
        length=p.size()
        sfxn = ScoreFunctionFactory.create_score_function('beta')
        sfxn(p)
        start_position=40
        end_position=length-40
        nc_high=0
        nc_low=999
        for ires in range(start_position,end_position,5):
            nc=count_contacts_accross_junction(p,ires)
            if(nc<nc_low and nc>0):
                nc_low=nc
            if(nc>nc_high):
                nc_high=nc
        output.write("{},{},{}\n".format(nc_low,nc_high,pdb_name))
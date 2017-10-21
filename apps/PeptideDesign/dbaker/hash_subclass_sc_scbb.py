#!/software/miniconda3/envs/pyrosetta3/bin/python
import pyrosetta
from pyrosetta import rosetta
from collections import defaultdict
#from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from rosetta.numeric import xyzVector_double_t as V3
from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
#from rif.hash import RosettaStubHash
import rif.hash
from rosetta import *
import pickle
import math
from collections import namedtuple


#  bb_hash has basic hash functions and bin resolution parameters

#  make_hash is used in generating hash from MC

#  use_hash is used in identifying hits in provided poses



def generate_canonical_rotamer_residues( residue_name , target_phi_psi):
    canonical_phi_psi = {"helical" : (-66.0,-42.0), "sheet" : (-126.0,124.0)}
    test_sequence = "AAX[%s]AA" % residue_name
    if target_phi_psi in canonical_phi_psi:
        target_phi , target_psi = canonical_phi_psi[target_phi_psi]
        print(target_phi)
        print(target_psi)
    else:
        target_phi,target_psi = target_phi_psi
    sf = pyrosetta.get_score_function()
    tryrot = TryRotamers(3, sf, 0, 0, True, solo_res=False, include_current=False )
    test_pose = pyrosetta.rosetta.core.pose.Pose()
    pyrosetta.rosetta.core.pose.make_pose_from_sequence( test_pose, test_sequence, "fa_standard" )
    for i in range(1,test_pose.size()+1):
        test_pose.set_psi(i, target_psi)
        test_pose.set_phi(i, target_phi)
        test_pose.set_omega(i, 180)
    
    tryrot.setup_rotamer_set( test_pose )
    rotamer_set = tryrot.rotamer_set()
    rotamers = [rotamer_set.rotamer(i).clone() for i in range(1, rotamer_set.num_rotamers() + 1)]
    return rotamers


# hash is from sc of sc_res to sc and bb of scbb_res.   key depends only on transform and possibly phi and psi of scbb_res.    Hence only require frame for res1, but pose for res2
# info stored may be res_names, chis of both residues, and relevant phipsi of scbb_res 
class bb_hash_sc_scbb:
    def __init__(self,hash_function,scbb_torsions):
        self.hash_function= hash_function
        self.scbb_torsions = scbb_torsions
#        self.store_res_names = store_res_names
#        self.name_data=name_data

    def get_bin_from_frame_and_pose(self, frame, pose, resN):
        scbb_res=pose.residue(resN)
        Aatom=frame.global2local(scbb_res.xyz('N'))
        Batom=frame.global2local(scbb_res.xyz('CA'))
        Catom=frame.global2local(scbb_res.xyz('C'))
        bb_stub = rosetta.core.kinematics.Stub(Aatom, Batom, Catom)
        tors=[]
# the following should work generally for the 0, 1 and 2 torsion dependency cases
        for x in self.scbb_torsions:
            if x in ['phi','psi']: tors.append( getattr(pose(resN),x))   # don't hash on chi1 as won't know it when scanning pose for hits
        k=self.hash_function(bb_stub,*tors)  # need to spell out items in list
        return k
            
    def get_frame_from_res(self,res):
        prot_bb=V3(res.xyz('N')),V3(res.xyz('CA')),V3(res.xyz('C'))
        frame = rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
        return frame

# we are hashing sc to scbb interactions (mainly hbonds).  will keep sc_res fixed, and sample rigid body orientation and relevant backbone torsions of scbb_res.
# so just have to orient rotamers once onto sc res, and then keep them fixed. 
class make_hash_sc_scbb(bb_hash_sc_scbb):
    def __init__(self,hash_function, scbb_torsions, orient_atoms, inverse_rot_backbones, sc_res, store_res_names_in_hash,name_data=[] ):
        bb_hash_sc_scbb.__init__(self,hash_function,scbb_torsions)
        if store_res_names_in hash:
            self.name_data=name_data
        self.n_added=0
        v=pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t()
        for atom in orient_atoms:
            v.append( (atom, atom) )
            
        self.dd=defaultdict(set)
            
 # generate and store all inverse rotamers of sc_res, their frames, and their chis            
        rots=get_rots(inverse_rot_backbones) 
        frames=[]
        chis=[]
        rot_info=[]
        for rot in rots:
            rot.orient_onto_residue(sc_res,self.v)
            frame=self.get_frame_from_res(rot)
            chi_list=[]
            for i in rot.nchi()+1:
                chi_list.append(int(rot.chi(i))/10)
            chi_tuple=tuple(chi_list)
            rot_data={"rot" : rot, "frame" : frame,"chis" : chi_tuple}
            rot_info.append(rot_data)
        self.rot_list=rot_info


    def get_rots(inverse_rot_backbones):
         rots=[]
         for angles in inverse_rot_backbones:
            rots = rots + generate_canonical_rotamer_residues(self.sc_res.name3(),angles)
         return rots
        
    def update_hash(self,pose,scbb_resN):
        store_data = namedtuple('store_data', ['res_names','sc_res_chis','scbb_res_chis']) # don't need to store scbb_tors as these are specified in key
        for rot_info in self.rot_list: 
            k=self.get_bin_from_frame_and_pose_bb(rot_info["frame"],pose,scbb_resN)
            if self.store_res_names : store_data.res_names=self.name_data
            store_data.sc_res_chis=rot_info['chis']
            scbb_res_chis=[]
            if chi1 in self.scbb_torsions: scbb_res_chis.append(int(pose.residue(scbb_resN).chi(1))/10)
            if chi2 in self.scbb_torsions: scbb_res_chis.append(int(pose.residue(scbb_resN).chi(2))/10)
            store_data.scbb_res_chis=tuple(scbb_res_chis)
            self.dd[k].add( store_data )
            self.n_added+=1
            if self.n_added % 1000==0:
                print('added %s hash size %s'%(self.n_added,len(self.dd.keys())))

#not used currently, use if only want to store most common rot
    def find_most_frequent_rot(self):
        d1={}
        for bin in dd.keys():
            v=list(dd[bin].values())
            k=list(dd[bin].keys())
            common_rot= k[v.index(max(v))]
            print(bin,k,v,common_rot)
        d1[bin]=common_rot
        return d1

class use_hash_sc_scbb(bb_hash_sc_scbb):
    def __init__(self,hash_function,scbb_torsions,filename,res_names_in_hash,name_data=[]):
        bb_hash_sc_scbb.__init__(self,hash_function,scbb_torsions)
        if not res_names_in hash:self.name_data=name_data  # if names are not in hash, need to specify them here
        self.dd =pickle.load(open(filename,"rb"))
        self.observed_keys=list(self.dd.keys())
        
 #two options when using hash:  1) count number of hits (bidentate hbs) and replace residue types in pose for each hit and
# 2) return full info for each hit (all rots and relevant chis).            
    def count_sc_scbb(self,pept,prot, contact_list):
#        mr = protocols.simple_moves.MutateResidue()
        nhb=0
        res_pair=namedtuple('res_pair',['pept_resN','prot_resN','hash_info'])
        res_pair_list=[]

        for resN in contact_list:
             frame = self.get_frame_from_res(prot.residue(resN))
             for i in range(2,pept.size()):
                k=self.get_bin_from_frame_and_pose(self, frame, pept, i)
                if k in self.observed_keys:
                  nhb=nhb+1
                  hash_info=next(iter(self.dd[k]))   # here just take random rot in hash bin
                  res_pair_list.append( res_pair(i,resN,hash_info))

        if nhb > 0: self.replace_residues_and_chis(pept,prot,res_pair_list)
        return nhb
    
#use resN for number, resName for name to avoid confusion
    def replace_residues_and_chis(self,pept,prot,res_pair_list):
         mr = protocols.simple_moves.MutateResidue()
         for res_pair in res_pair_list:
             #need to unpack residue names, the relevant backbone torsions, and the chi angles
             # residue names are either stored in hash, or specified when hash object is set up. advantage of former is that same hash can contain multiple different res pairs, the advantage of latter is that it is more memory efficient.
             if self.res_names_in_hash:
                 mr.set_res_name(res_pair.hash_info.sc_resName)
             else:
                 mr.set_res_name(self.name_data.sc_resName)
             mr.set_target(res_pair.prot_resN)
             mr.apply(prot)
             nchi=1  #chis start at 1
             for chi in res_pair.hash_info.sc_res_chis():
                 prot.set_chi(nchi,res_pair.prot_resN,chi)
                 nchi=nchi+1
             if self.res_names_in hash:
                 mr.set_res_name(res_pair.hash_info,scbb_resName)
             else:
                 mr.set_res_name(self.name_data.scbb_resName)
             mr.set_target(res_pair.pept_resN)
             mr.apply(prot)
             nchi=1  #chis start at 1
             for chi in res.hash_info.scbb_res_chis():
                 prot.set_chi(nchi,res_pair.prot_pept_resN,chi)
                 nchi=nchi+1

    def get_all_rots(self,pept,prot,contact_list):
        pept_rots={}
        prot_rots={}
        nhb=0
        for resN in contact_list:
             frame = self.get_frame_from_res(prot.residue(resN))
             for i in range(2,pept.size()):
                k=self.get_bin_from_frame_and_pose(self, frame, pept, i)
                if k in list(self.dd.keys()):
                  nhb=nhb+1
                  for store_data in self.dd[k]:
                      prot_rots[resN].append(store_data.sc_res_chis())
                      pept_rots[i].append(store_data.scbb_res_chis())
        return nhb,pept_rots,prot_rots
                 
    def convert_to_set(self):
        key_set=set(self.dd.keys())
        return key_set



    # def count_nhb(self,pept,prot):
    #     nhb=0
    #     prot_res=prot.residue(1) #skip termini as atom types differ
    #     pept_residue=pept.residue(1)
    #         #for prot_res in prot_residues:
    #     frame = self.get_frame_from_res(prot_res)
    #     bb_rays=get_bb_rays_from_res(pept_residue)
    #     k=self.get_bin_from_rot_and_frame(bb_rays,frame)
    #     if k in self.dd.keys():
    #         nhb=nhb+1
    #         chis=next(iter(self.dd[k]))   # for now, use random rot
    #         self.res_list.append( (prot_res.seqpos(),10*chis[0],10*chis[1]) )
    #         print prot_res.seqpos(),prot_res.name1(),pept_residue.seqpos(),pept_residue.name1()
    #     #if nhb>0: self.replace_residues_and_chis(prot)
        
    #     return nhb
    

    
    # def replace_residues_and_chis(self,pept,prot):
    #      mr = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
    #      mr.set_res_name('HIS')
    #      for res_info in self.res_list:
    #          mr.set_target(res_info[0])
    #          mr.apply(prot)
    #          prot.set_chi(1,res_info[0],res_info[1])
    #          prot.set_chi(2,res_info[0],res_info[2])
             


    # def get_all_rots(self,pept,prot):

    #     rots={}
    #     nhb=0
    #     prot_residues=[prot.residue(i) for i in range(2,prot.size())] #skip termini as atom types differ
    #     pept_residues=pept.residue(1)
    #     self.res_list=[]
    #     for prot_res in prot_residues:
    #          frame = self.get_frame_from_res(prot_res)
    #          for pept_residue in pept_residues:
    #             bb_rays=self.get_bb_rays_from_res(pept_residue)
    #             k=self.get_bin_from_rot_and_frame(bb_rays,frame)
    #             if k in self.dd.keys():
    #               nhb=nhb+1
    #               rots[prot_res.seqpos()]=self.dd[k]
    #               print prot_res.seqpos(),prot_res.name1(),pept_residue.seqpos(),pept_residue.name1()
        
    #     return nhb,rots

    def test_hash(self,pdb,sc_resN,scbb_resN):
        p=rosetta.core.import_pose.pose_from_file(pdb)
        bb_res=p.residue(2)
        asn_res=p.residue(4)
        frame = self.get_frame_from_res(pose.residue(sc_resN))
        k=self.get_bin_from_frame_and_pose(self, frame, p, sc_resN)                                        
        bb_rays=self.get_bb_rays_from_res(pept_residue)
        if k in self.dd.keys():
            print(dd[k])
            return 1
        else:
            return 0

            
import string
if __name__ == "__main__":
    pyrosetta.init()
    prot=rosetta.core.import_pose.pose_from_file("5AEI_A.pdb")
    pept=rosetta.core.import_pose.pose_from_file("5AEI_D.pdb")
#    hash_chi=use_hash_rot(1.0,3.0,"bb_asn_hash_combo_new_ray_1.0_rot0")
#    hash_nc=use_hash(1.0,3.0,"bb_asn_hash_combo_new_ray_1.0_rot0")
    hash_nc=use_hash_rot(1.0,3.0,"bb_asn_dict_combo_1.0_rot0")
    nhb=hash_nc.count_asn_bb(pept,prot)
#    s=hash_nc.convert_to_set()
#    pickle.dump(s,open("bb_asn_dict_combo_1.0_set","wb"))
    print(nhb)
    pdb_list=map(string.split,open('pdb.list','r').readlines())
    for pdb in pdb_list:
        p=rosetta.core.import_pose.pose_from_file(pdb[0])
        pept,prot=p.split_by_chain()
        nhb=hash_nc.count_asn_bb(pept,prot)
        print(nhb,pdb)
        


#!/software/miniconda3/envs/pyrosetta3/bin/python
import pyrosetta
from pyrosetta import rosetta
from collections import defaultdict
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from rosetta.numeric import xyzVector_double_t as V3
from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
from rif.hash import RosettaStubHash
from rosetta import *
import pickle
import math

#  bb_hash has basic hash functions and bin resolution parameters

#  make_hash adds in rotamer superposition
#  make_hash_no_rot  and make_hash_store_rot  add in table updating wih set and defaultdict respectively.  make_hash_store_and_count_rot keeps only most frequent rotamer

#  use_hash inherits from bb_hash--identifies hits in pdb
#  use_hash_rot  same, but copies in rot info also


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


# assume that hash will depend only on transform and phi and psi of second residue.  first residue is a sidechain, second is either or both sidechain and backbone (so might need phi and psi).  Hence only require frame for res1, but pose for res2

class bb_hash_6-8D:
    def __init__(self,resl,lever,hash_function,res2_key_vars,res1_store_vars,res2_store_vars):
        self.resl=resl
        self.lever=lever
        self.binner= hash_function(resl, lever)
        self.res2_key_vars = res2_key_vars
        self.res1_store_vars = res1_store_vars
        self.res2_store_vars = res2_store_vars

    def get_bin_from_frame_and_pose(self, frame, pose, resN):
#        prot_bb=V3(res1.xyz('N')),V3(res1.xyz('CA')),V3(res1.xyz('C'))
#        frame = rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
        
        res2=pose.residue(resN)
        Aatom=frame.global2local(res2.xyz('N'))
        Batom=frame.global2local(res2.xyz('CA'))
        Catom=frame.global2local(res2.xyz('C'))
        bb_stub = rosetta.core.kinematics.Stub(Aatom, Batom, Catom)
        tors=[]
        for x in self.res2_key_vars:
            tors.append( getattr(pose(resN),x))
        k=self.hash_function(bb_stub,*tors)  # need to spell out items in list
        return k
        
    
    def get_frame_from_res(self,res):
        prot_bb=V3(res.xyz('N')),V3(res.xyz('CA')),V3(res.xyz('C'))
        frame = rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
        return frame

class make_hash(bb_hash_6-8D):
    def __init__(self,resl,lever,hash_function,res2_key_vars,res1_store_vars,res2_store_vars, orient_atoms, base_res):
        bb_hash_6D-8D.__init__(self,resl,lever,hash_function,res2_key_vars,res1_store_vars,res2_store_vars)
        self.n_added=0
        self.rots = []
        self.start_frame=get_frame_from_res(base_res)  # residue that is fixed in MC
        v=pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t()
        for atom in orient_atoms:
            v.append( (atom, atom) )
        self.v=v

        #self.rots= generate_canonical_rotamer_residues("HIS","helical")
        
    def orient_rots(self,base,name):
        new_rotamers = []
        self.rots = []
        for angles in self.target_phi_psi:
            self.rots = self.rots + generate_canonical_rotamer_residues(self.res1.name3(),angles)
        for rot in self.rots:
            rot.orient_onto_residue(base,self.v)
            new_rotamers.append(rot)
        return new_rotamers
'''
class make_hash_no_rot(make_hash):
    def __init(self,resl,lever):
        make_hash.__init__(self,resl,lever)
        self.s=set()
        
    def update_table(self,gn):
        orient_rots(gn.residue(4))  # rots is transformed
        bb=gn.residue(2)
        bb_rays=get_bb_rays_from_res(bb)
        for rot in self.rots:
            k=self.get_bin_from_rot_and_bb(self,bb_rays,rot)
            self.s.add(k)
            self.n_added+=1
            if self.n_added % 1000==0:
                print('added %s hash size %s'%(self.n_added,len(s)))
'''
class make_hash_store_rot(make_hash):
    def __init__(self,resl,lever,target_angle,orient_atoms):
        make_hash.__init__(self,resl,lever,target_angle,orient_atoms)
        self.dd=defaultdict(set)
        self.point = int(-1000)
        self.point2 = int(0)
        self.mark = 'False'
        
    def update_table_chi(self,bb_rays,sample_resi,resi_name,chi_number):
        # orient rotamers onto the residue conformation which pass through MC
        new_rots = self.orient_rots(sample_resi,resi_name)  # rots is transformed
        #bb=gn.residue(1)
        for rot in new_rots:
            chi_value=[]
            for i in range(1,int(chi_number+1)):
                chi_value.append(int(rot.chi(i)/10))
            #chi1=int(rot.chi(1)/10)
            #print chi1
            #chi2=int(rot.chi(2)/10)
            #print chi2
            k=self.get_bin_from_rot_and_bb(bb_rays,rot)
            #print len(tuple(chi_value))
            self.dd[k].add( tuple(chi_value) )
            self.n_added+=1
            if self.n_added % 1000==0:
                print('added %s hash size %s'%(self.n_added,len(self.dd.keys())))
                self.point = self.point2
                self.point2 = len(self.dd.keys())
                if self.point2 - self.point < 5:
                    self.mark = 'True'
        return self.mark
'''
class make_hash_store_and_count_rot(make_hash):
    def __init__(self,resl,lever):
        make_hash.__init__(self,resl,lever)
        self.dd=defaultdict(lambda : defaultdict(int))
        
    def update_table_chi(self,gn,rots):
        orient_rots(gn.residue(4))  # rots is transformed
        bb=gn.residue(2)
        bb_rays=self.get_bb_rays_from_res(bb)
        for rot in rots:
            chi1=int(rot.chi(1)/10)
            chi2=int(rot.chi(2)/10)
            k=self.get_bin_from_rot_and_bb(self,bb_rays,rot)
            self.dd[k][(chi1,chi2)] += 1
            self.n_added+=1
            if self.n_added % 1000==0:
              print('added %s hash size %s'%(n_added,len(dd.keys())))

    def find_most_frequent_rot(self):
        d1={}
        for bin in dd.keys():
            v=list(dd[bin].values())
            k=list(dd[bin].keys())
            common_rot= k[v.index(max(v))]
            print bin,k,v,common_rot
        d1[bin]=common_rot
        return d1

class use_hash(bb_hash):
    def __init__(self,resl,lever,filename):
        bb_hash.__init__(self,resl,lever)
        self.s =pickle.load(open(filename,"rb"))

    def replace_residues(self,pose):
         mr = protocols.simple_moves.MutateResidue()
         mr.set_res_name('ASN')
         for res in self.res_list:
             mr.set_target(res)
             mr.apply(pose)
             
    def count_asn_bb(self,pept,prot):
        mr = protocols.simple_moves.MutateResidue()
        nhb=0
        prot_residues=[prot.residue(i) for i in range(2,prot.size())] #skip termini as atom types differ
        pept_residues=[pept.residue(i) for i in range(2,pept.size())]
        self.res_list=[]
        for prot_res in prot_residues:
             frame = self.get_frame_from_res(prot_res)
             for pept_residue in pept_residues:
                bb_rays=self.get_bb_rays_from_res(pept_residue)
                k=self.get_bin_from_rot_and_frame(bb_rays,frame)
                if k in self.s:
                  nhb=nhb+1
                  print prot_res.seqpos(),prot_res.name1(),pept_residue.seqpos(),pept_residue.name1()
                  self.res_list.append(prot_res.seqpos())
        if nhb>0: self.replace_residues(prot)
        
#                  mr.set_target(resn)
#                  mr.set_res_name('ASN')
#                  mr.apply(prot)
        return nhb


class use_hash_rot(bb_hash):
    def __init__(self,resl,lever,filename):
        bb_hash.__init__(self,resl,lever)
        self.dd =pickle.load(open(filename,"rb"))
        self.res_list=[]
        tt=self.dd.keys()
        #for key in tt[:5]:
            #print key, self.dd[key]
            
    def convert_to_set(self):
        key_set=set(self.dd.keys())
        return key_set
            
    def replace_residues_and_chis(self,prot):
         mr = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
         mr.set_res_name('HIS')
         for res_info in self.res_list:
             mr.set_target(res_info[0])
             mr.apply(prot)
             prot.set_chi(1,res_info[0],res_info[1])
             prot.set_chi(2,res_info[0],res_info[2])
             
    def count_protein_pept(self,pept,prot):
        nhb=0
        prot_res=prot.residue(1) #skip termini as atom types differ
        pept_residue=pept.residue(1)
            #for prot_res in prot_residues:
        frame = self.get_frame_from_res(prot_res)
        bb_rays=get_bb_rays_from_res(pept_residue)
        k=self.get_bin_from_rot_and_frame(bb_rays,frame)
        if k in self.dd.keys():
            nhb=nhb+1
            chis=next(iter(self.dd[k]))   # for now, use random rot
            self.res_list.append( (prot_res.seqpos(),10*chis[0],10*chis[1]) )
            print prot_res.seqpos(),prot_res.name1(),pept_residue.seqpos(),pept_residue.name1()
        #if nhb>0: self.replace_residues_and_chis(prot)
        
        return nhb

    def get_all_rots(self,pept,prot):
        rots={}
        nhb=0
        prot_residues=[prot.residue(i) for i in range(2,prot.size())] #skip termini as atom types differ
        pept_residues=pept.residue(1)
        self.res_list=[]
        for prot_res in prot_residues:
             frame = self.get_frame_from_res(prot_res)
             for pept_residue in pept_residues:
                bb_rays=self.get_bb_rays_from_res(pept_residue)
                k=self.get_bin_from_rot_and_frame(bb_rays,frame)
                if k in self.dd.keys():
                  nhb=nhb+1
                  rots[prot_res.seqpos()]=self.dd[k]
                  print prot_res.seqpos(),prot_res.name1(),pept_residue.seqpos(),pept_residue.name1()
        
        return nhb,rots

    def test_hash(self,pdb):
        p=rosetta.core.import_pose.pose_from_file(pdb)
        bb_res=p.residue(2)
        asn_res=p.residue(4)
        frame = self.get_frame_from_res(asn_res)
        bb_rays=self.get_bb_rays_from_res(pept_residue)
        k=self.get_bin_from_rot_and_frame(bb_rays,frame)
        if k in self.dd.keys():
            print dd[k]
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
    print nhb
    pdb_list=map(string.split,open('pdb.list','r').readlines())
    for pdb in pdb_list:
        p=rosetta.core.import_pose.pose_from_file(pdb[0])
        pept,prot=p.split_by_chain()
        nhb=hash_nc.count_asn_bb(pept,prot)
        print nhb,pdb
        
    #    prot.dump_pdb('%s.pdb'%nhb)
    

'''

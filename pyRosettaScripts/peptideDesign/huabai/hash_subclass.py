#!/home/sheffler/venv/david/bin/ipython
import pyrosetta
from pyrosetta import rosetta
from collections import defaultdict
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from rosetta.numeric import xyzVector_double_t as V3
from rosetta import *
import pickle
import math

#  bb_hash has basic hash functions and bin resolution parameters

#  make_hash adds in rotamer superposition
#  make_hash_no_rot  and make_hash_store_rot  add in table updating wih set and defaultdict respectively.  make_hash_store_and_count_rot keeps only most frequent rotamer

#  use_hash inherits from bb_hash--identifies hits in pdb
#  use_hash_rot  same, but copies in rot info also

class bb_hash:
    def __init__(self,resl,lever):
        self.resl=resl
        self.lever=lever
        self.binner= RayRay10dHash(resl, lever)

    def get_frame_from_res(self,res):
        prot_bb=V3(res.xyz('N')),V3(res.xyz('CA')),V3(res.xyz('C'))
        frame = rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
        return frame

    def get_bb_rays_from_res(self,res):
        bb_rays= V3(res.xyz('C')),V3(res.xyz('O')),V3(res.xyz('N')),V3(res.xyz('H'))
        # bb_rays= V3(res.xyz('C')),V3(res.xyz('O')),V3(res.xyz('CA')),V3(res.xyz('N'))
        return bb_rays

    def get_bin_from_rot_and_bb(self,bb_rays,rot):
        frame=get_frame_from_res(rot)
        Catom=frame.global2local(bb_rays[0])
        Oatom=frame.global2local(bb_rays[1])
        Natom=frame.global2local(bb_rays[2])
        Hatom=frame.global2local(bb_rays[3])
        COray = Ray(orig=Oatom, dirn=(Oatom-Catom))
        NHray = Ray(orig=Hatom, dirn=(Hatom-Natom))
        k=self.binner.get_key(COray,NHray)
        return k

    def get_bin_from_rot_and_frame(self,bb_rays,frame):    #use this one when have many chain 2 residues for each chain 1 residue; can precompute frame from chain 1 res

        Catom=frame.global2local(bb_rays[0])
        Oatom=frame.global2local(bb_rays[1])
        Natom=frame.global2local(bb_rays[2])
        Hatom=frame.global2local(bb_rays[3])
        COray = Ray(orig=Oatom, dirn=(Oatom-Catom))
        NHray = Ray(orig=Hatom, dirn=(Hatom-Natom))
        k=self.binner.get_key(COray,NHray)
        return k

class make_hash(bb_hash):
    def __init__(self,resl,lever):
        bb_hash.__init__(self,resl,lever)
        self.n_added=0
        self.rots= generate_canonical_rotamer_residues("ASN","helical")
        v=pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t()
        v.append( ('OD1', 'OD1') )
        v.append( ('ND2', 'ND2') )
        v.append( ('CG', 'CG') )
        self.v=v

    def orient_rots(base):
        for rot in self.rots:
            rot.orient_onto_residue(base,self.v)

class make_hash_no_rot(make_hash):
    def __init(self,resl,lever):
        make_hash.__init__(self,resl,lever)
        self.s=set()

    def update_table(self,gn):
        orient_rots(gn.residue(4))  # rots is transformed
        bb=gn.residue(2)
        bb_rays=self.get_bb_rays_from_res(bb)
        for rot in self.rots:
            k=self.get_bin_from_rot_and_bb(self,bb_rays,rot)
            self.s.add(k)
            self.n_added+=1
            if self.n_added % 1000==0:
                print('added %s hash size %s'%(self.n_added,len(s)))

class make_hash_store_rot(make_hash):
    def __init__(self,resl,lever):
        make_hash.__init__(self,resl,lever)
        self.dd=collections.defaultdict(set)

    def update_table_chi(self,gn,rots):
        orient_rots(gn.residue(4))  # rots is transformed
        bb=gn.residue(2)
        bb_rays=self.get_bb_rays_from_res(bb)
        for rot in rots:
            chi1=int(rot.chi(1)/10)
            chi2=int(rot.chi(2)/10)
            k=self.get_bin_from_rot_and_bb(self,bb_rays,rot)
            self.dd[k].add( (chi1,chi2) )
            self.n_added+=1
            if self.n_added % 1000==0:
              print('added %s hash size %s'%(n_added,len(dd.keys())))

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
         # mr.set_res_name('ASN')
         mr.set_res_name('TYR')
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
        tt=self.dd.keys()
        for key in tt[:5]:
            print key, self.dd[key]

    def convert_to_set(self):
        key_set=set(self.dd.keys())
        return key_set

    def replace_residues_and_chis(self,prot):
         mr = protocols.simple_moves.MutateResidue()
         mr.set_res_name('TYR')
         for res_info in self.res_list:
             mr.set_target(res_info[0])
             mr.apply(prot)
             prot.set_chi(1,res_info[0],res_info[1])
             prot.set_chi(2,res_info[0],res_info[2])

    def count_asn_bb(self,pept,prot):
        nhb=0
        prot_residues=[prot.residue(i) for i in range(2,prot.size())] #skip termini as atom types differ
        pept_residues=[pept.residue(i) for i in range(2,pept.size())]
        self.res_list=[]
        for prot_res in prot_residues:
             frame = self.get_frame_from_res(prot_res)
             for pept_residue in pept_residues:
                bb_rays=self.get_bb_rays_from_res(pept_residue)
                k=self.get_bin_from_rot_and_frame(bb_rays,frame)
                if k in self.dd.keys():
                  nhb=nhb+1
                  chis=next(iter(self.dd[k]))   # for now, use random rot
                  self.res_list.append( (prot_res.seqpos(),10*chis[0],10*chis[1]) )
                  print prot_res.seqpos(),prot_res.name1(),pept_residue.seqpos(),pept_residue.name1()
        if nhb>0: self.replace_residues_and_chis(prot)

        return nhb

    def get_all_rots(self,pept,prot):
        rots={}
        nhb=0
        prot_residues=[prot.residue(i) for i in range(2,prot.size())] #skip termini as atom types differ
        pept_residues=[pept.residue(i) for i in range(2,pept.size())]
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


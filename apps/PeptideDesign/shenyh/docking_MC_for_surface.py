#!/home/sheffler/venv/david/bin/ipython
import pyrosetta
from xyzMath import *
from extract_silica_unit import *
from hash_subclass import *
import random
import time
import math
def center_on_z_axis(res1,res2,pose):
    stub1 = stub(res1[0],res1[1],res1[2])
    stub2 = stub(res2[0],res2[1],res2[2])
    xform = stub2 * ~stub1
    axis, ang, cen = xform.rotation_axis_center()
    rad_vec= projperp(axis, res1[1] - cen).normalized()
    yvec=axis.cross(rad_vec)
    stub1=stub(cen+rad_vec,cen+yvec,cen+axis)
    stub2=stub( Vec(1,0,0),Vec(0,1,0),Vec(0,0,1) )
    transform=stub2*~stub1
    #    transform= Xform().from_two_vecs(axis,cen-res1[0])
    
    for ir in range(1, pose.size() + 1):
        for ia in range(1, pose.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(pose.xyz(aid))
            newxyz = transform * oldxyz
            pose.set_xyz(aid, newxyz.to_rosetta())
    return center_com_z(pose.clone())

def center_com_z(p):
    z_tot=0.
    n=0
    for ir in range(1,p.size()+1):
        for ia in range(1,p.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            z = p.xyz(aid).z
            n=n+1
            z_tot=z_tot+z
    z_vec=Vec(0.,0.,z_tot/float(n))
    #   print('z offset %s'%z_vec)
    for ir in range(1, p.size() + 1):
        for ia in range(1, p.residue(ir).natoms() + 1):
            aid = pyrosetta.rosetta.core.id.AtomID(ia, ir)
            newxyz = Vec(p.xyz(aid))-z_vec
            p.set_xyz(aid, newxyz.to_rosetta())
    return p

pyrosetta.init("-extra_res_fa OSi.params Si.params")
p=pyrosetta.rosetta.core.import_pose.pose_from_file("final_scaffold_onsurface.pdb")
sf=pyrosetta.get_score_function()
chains=p.split_by_chain()
scaffold=chains[3]
old_E=sf(p)
low_E=old_E
low_E_num=0
mute_list = [x for x in range(1,scaffold.size()+1) if x%18 == 13]
mute = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
mute.set_res_name('ALA')
for i in mute_list:
    mute.set_target(i)
    mute.apply(scaffold)

scaffold.dump_pdb("docking_origin.pdb")
surface=chains[1]
surface_1=chains[2]
surface.append_pose_by_jump(surface_1,50)
complex=surface.clone()
complex.append_pose_by_jump(scaffold.clone(),100)
complex.dump_pdb("complex.pdb")
surface_units = get_surface_unit(surface)
prot=complex.jump(100)
total_cycle = int(100)
T=0.2
n_accept=0
print "origin_energy={}".format(old_E)
print "read hashtable"
t0=time.time()
hash_nc = use_hash_rot(0.4,3.0,"SiO2_hashtable_all")
nsurf = hash_nc.count_protein_pept(scaffold,surface,surface_units)
print nsurf
print time.time()-t0
mm = pyrosetta.rosetta.core.kinematics.MoveMap()
mm.set_chi_true_range(101,244)
if nsurf>0:
    scaffold.dump_pdb("origin_scaffold_design.pdb")
print "begin MC"
for i in range(total_cycle):
    print time.time()-t0
    rand=random.random()
    k=1
    if rand <0.5 :
        k = -1
    prot=complex.jump(100)
    prot.gaussian_move(k,0.15,0.5)
    complex.set_jump(100,prot)
    #complex.dump_pdb("test_gas_{}.pdb".format(i))
    chains=complex.clone().split_by_chain()
    scaffold=chains[3]
    for j in mute_list:
        mute.set_target(j)
        mute.apply(scaffold)
    nsurf = hash_nc.count_protein_pept(scaffold,surface,surface_units)
    print nsurf
    if nsurf>0:
        scaffold.dump_pdb("{}_scaffold_design.pdb".format(i))
        complex_new=surface.clone()
        complex_new.append_pose_by_jump(scaffold.clone(),100)
        min_mover = pyrosetta.rosetta.protocols.simple_moves.MinMover(mm,sf,'dfpmin_armijo_nonmonotone',0.01,True)
        min_mover.apply(complex_new)
        new_E=sf(complex_new)
        print new_E
        delta = new_E-old_E
        accept=1
        if delta>0:
            if math.exp(-delta/T)<random.random():
                accept = 0
        if accept:
            old_E = new_E
            complex = complex_new
            if low_E>new_E:
                low_E = new_E
                low_E_num = i
            if new_E < 2500:
                #n_accept = n_accept + 1
                complex_new.dump_pdb("output_complex_{}.pdb".format(i))

print "low_energy={}".format(low_E)
print "low_energy_conformation=output_complex_{}.pdb".format(low_E_num)


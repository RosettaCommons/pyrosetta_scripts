#!/home/sheffler/venv/david/bin/ipython
from xyzMath import *
from hash_subclass import *
from pyrosetta import *
from extract_silica_unit import *
import random

def rotate_around_axis(axis, angle, pose):
    r_axis=axis.to_rosetta()
    t=rotation_matrix_degrees(axis,angle)
    for ir in range(1,pose.size() + 1):
        for ia in range(1,pose.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia,ir)
            oldxyz = Vec(pose.xyz(aid))
            newxyz = t * oldxyz
            pose.set_xyz(aid, newxyz.to_rosetta())
    return pose

def translate_along_axis(dis_vec,pose):
    for ir in range(1,pose.size() + 1):
        for ia in range(1,pose.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia,ir)
            newxyz = Vec(pose.xyz(aid)) + dis_vec
            pose.set_xyz(aid,newxyz.to_rosetta())
    return pose


def dock_peptide(p,nbins):
    for x in range(nbins):
        for y in range(nbins):
            for z in range(nbins):
                disx = -0.5 + 1.0 * float(x)/float(nbins)
                disy = -1.0 + 1.0 * float(y)/float(nbins)
                disz = -0.5 + 1.0 * float(z)/float(nbins)
                b4 = translate_along_axis( Vec(float(disx),0.0,0.0),p.clone())
                b5 = translate_along_axis( Vec(0.0,float(disy),0.0),b4.clone())
                b6 = translate_along_axis( Vec(0.0,0.0,float(disz)),b5.clone())
                yield(b6)

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
prot = pyrosetta.rosetta.core.import_pose.pose_from_file("mute_6k+1_6k+5.pdb")
surface = pyrosetta.rosetta.core.import_pose.pose_from_file("silica_surface_near_protein.pdb")
'''
mute = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
mute.set_res_name('ALA')
for i in range(1,prot.size()+1):
    if i % 6 ==1 or i % 6 ==5 :
        mute.set_target(i)
        mute.apply(prot)
prot.dump_pdb("mute_6k+1_6k+5.pdb")

res1=[Vec(prot.residue(8).xyz("N")),Vec(prot.residue(8).xyz("CA")),Vec(prot.residue(8).xyz("C"))]
res2=[Vec(prot.residue(14).xyz("N")),Vec(prot.residue(14).xyz("CA")),Vec(prot.residue(14).xyz("C"))]
prot=center_on_z_axis(res1,res2,prot)
hash_nc = use_hash_rot(0.4,3.0,"CaCO3_hashtable")
#prot = pyrosetta.rosetta.core.import_pose.pose_from_file("hisonly1.pdb")
#pept = pyrosetta.rosetta.core.import_pose.pose_from_file("hemeonly1.pdb")
#pept.dump_pdb("testifheme.pdb")
#prot.dump_pdb("testifhis.pdb")
'''
surface_units = get_surface_unit(surface)
print surface_units

hash_nc = use_hash_rot(0.4,3.0,"SiO2_onelys_allsurface_hashtable")
#num = 0
rotamer_gen = dock_peptide(prot,5)
for rotamer in rotamer_gen:
    #num = num +1
    #rotamer.dump_pdb("rotamer_test_{}.pdb".format(num))

    print 'Docking begins!'
    nsurf = hash_nc.count_protein_pept(rotamer,surface,surface_units)
    print nsurf
    #relax the whole structure and output the energy score to see if the structure can be accepted
    if nsurf > 0:
        rotamer.dump_pdb(("%s_his_protein.pdb" % nsurf))
    #relax the whole structure and output the energy score to see if the structure can be accepted



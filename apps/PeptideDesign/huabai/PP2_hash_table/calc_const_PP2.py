#!/home/sheffler/venv/david/bin/python
import pyrosetta
from xyzMath import *
from math import *
from rosetta.numeric import xyzVector_double_t as V3

def get_normal_vector(atom1, atom2, atom3):
    a1_x = atom1[0]
    a1_y = atom1[1]
    a1_z = atom1[2]
    a2_x = atom2[0]
    a2_y = atom2[1]
    a2_z = atom2[2]
    a3_x = atom3[0]
    a3_y = atom3[1]
    a3_z = atom3[2]
    a = ((a2_y - a1_y)*(a3_z - a1_z) - (a2_z - a1_z)*(a3_y - a1_y))
    b = ((a2_z - a1_z)*(a3_x - a1_x) - (a2_x - a1_x)*(a3_z - a1_z))
    c = ((a2_x - a1_x)*(a3_y - a1_y) - (a2_y - a1_y)*(a3_x - a1_x))
    return Vec(a,b,c)

def get_vector_angl(vec1,vec2):
    multi = v_normal * v_CH
    length = sqrt(float(v_normal * v_normal)) * sqrt(float(v_CH * v_CH))
    angl = acos(float(multi)/float(length))
    return angl

def rotation_around_axis(axis, angle, pose):
    r_axis=axis.to_rosetta()
    t=rotation_matrix_degrees(r_axis,angle)
    print t
    aid = pyrosetta.rosetta.core.id.AtomID(24,4)
    oldxyz = Vec(pose.xyz(aid))
    newxyz = t * oldxyz
    pose.set_xyz(aid, newxyz.to_rosetta())
    return pose
pyrosetta.init()
p = pyrosetta.rosetta.core.import_pose.pose_from_file("origin_complex.pdb")
q = p.clone()
prot_bb=V3(p.residue(4).xyz('CZ')),V3(p.residue(4).xyz('OH')),V3(p.residue(4).xyz('HH'))
frame = pyrosetta.rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
Hatom_local = Vec(frame.global2local(p.residue(4).xyz('HH')))
print Hatom_local
Catom_local = Vec(frame.global2local(p.residue(4).xyz('CZ')))
print Catom_local
Oatom_local = Vec(frame.global2local(p.residue(4).xyz('OH')))
print Oatom_local
v1 = Oatom_local
v2 = Catom_local
CO_vec = v1 - v2
print CO_vec
t=rotation_matrix_degrees(CO_vec,90.0)
Hatom_local_new = t*Hatom_local
Hatom_global_new = frame.local2global(Hatom_local_new.to_rosetta())
aid = pyrosetta.rosetta.core.id.AtomID(24,4)
q.set_xyz(aid, Hatom_global_new)
q.dump_pdb('rotate.pdb')




'''
pyrosetta.init()
p = pyrosetta.rosetta.core.import_pose.pose_from_file("origin_complex.pdb")
OH_line = p.xyz(pyrosetta.rosetta.core.id.AtomID(13,4)) - p.xyz(pyrosetta.rosetta.core.id.AtomID(12,4))
OH_vec = Vec(OH_line)
length = sqrt(float(OH_vec*OH_vec))
OH_vec = OH_vec/length
print OH_vec
p2 = rotation_around_axis(OH_vec, 90.0, p.clone())
p2.dump_pdb("rotate.pdb")
'''
'''
v_normal = get_normal_vector(p.xyz(pyrosetta.rosetta.core.id.AtomID(7,4)),p.xyz(pyrosetta.rosetta.core.id.AtomID(11,4)),p.xyz(pyrosetta.rosetta.core.id.AtomID(10,4)))
#1HDatom
v_H = Vec(p.xyz(pyrosetta.rosetta.core.id.AtomID(14,2))[0],p.xyz(pyrosetta.rosetta.core.id.AtomID(14,2))[1],p.xyz(pyrosetta.rosetta.core.id.AtomID(14,2))[2])
#CDatom
v_C = Vec(p.xyz(pyrosetta.rosetta.core.id.AtomID(7,2))[0],p.xyz(pyrosetta.rosetta.core.id.AtomID(7,2))[1],p.xyz(pyrosetta.rosetta.core.id.AtomID(7,2))[2])
v_CH = v_C - v_H
ang = get_vector_angl(v_normal , v_CH)
arom_center = p.xyz(pyrosetta.rosetta.core.id.AtomID(7,4)) + p.xyz(pyrosetta.rosetta.core.id.AtomID(11,4)) +p.xyz(pyrosetta.rosetta.core.id.AtomID(10,4))
arom_center_vec = Vec(float(arom_center[0])/3.0, float(arom_center[1])/3.0, float(arom_center[2])/3.0)
Hcenter_vec = v_H - arom_center_vec
Ccenter_vec = v_C - arom_center_vec
Hcenter_distance = sqrt(float(Hcenter_vec * Hcenter_vec))
Ccenter_distance = sqrt(float(Ccenter_vec * Ccenter_vec))
C_ang = get_vector_angl(v_normal , Ccenter_vec)
projection_distance = float(Ccenter_distance * sin(C_ang))
print projection_distance
print Hcenter_distance
print ang

                

multi = v_normal * v_CH
length = sqrt(v_normal * v_normal) * sqrt(v_CH * v_CH)
angl = acos(float(multi)/float(length))
print v


v1 = Vec(1.0,1.0,1.0)
v2 = Vec(2.0,2.0,2.0)
a = v1*v2
print a
'''

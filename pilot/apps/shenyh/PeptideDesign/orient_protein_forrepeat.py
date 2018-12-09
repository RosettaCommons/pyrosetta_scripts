#!/home/sheffler/venv/david/bin/ipython
import pyrosetta
from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
from xyzMath import *
import pyrosetta.toolbox.numpy_utils as np_utils
from numpy_utils import *
import time
import math

def get_vector_angl(vec1,vec2):
    multi = vec1 * vec2
    length = sqrt(float(vec1 * vec1)) * sqrt(float(vec2 * vec2))
    angl = acos(float(multi)/float(length))
    return angl

def generate_canonical_rotamer_residues( residue_name , target_phi_psi):
    canonical_phi_psi = {"helical" : (-66.0,-42.0), "sheet" : (-126.0,124.0)}
    test_sequence = "AAX[%s]AA" % residue_name
    if target_phi_psi in canonical_phi_psi:
        target_phi , target_psi = canonical_phi_psi[target_phi_psi]
        print target_phi
        print target_psi
    else:
        target_phi,target_psi = target_phi_psi
    sf = pyrosetta.get_score_function()
    tryrot = TryRotamers(3, sf, 2, 0, True, solo_res=False, include_current=False )
    test_pose = pyrosetta.rosetta.core.pose.Pose()
    pyrosetta.rosetta.core.pose.make_pose_from_sequence( test_pose, test_sequence, "fa_standard" )
    for i in range(1,test_pose.size()+1):
        test_pose.set_psi(i, target_psi)
        test_pose.set_phi(i, target_phi)
        test_pose.set_omega(i, 180)
    
    tryrot.setup_rotamer_set( test_pose )
    rotamer_set = tryrot.rotamer_set()
    rotamers = [rotamer_set.rotamer(i).clone() for i in xrange(1, rotamer_set.num_rotamers() + 1)]
    return rotamers

def get_anchor_coordinates_from_pose(residue):
    bb_atoms = ['N','CA','C']
    coords = []
    for atom in bb_atoms:
        coords.append([residue.xyz(atom).x,residue.xyz(atom).y,residue.xyz(atom).z])
    return np.mat(coords)

def get_vector_angl(vec1,vec2):
    multi = abs(vec1 * vec2)
    length = sqrt(float(vec1 * vec1)) * sqrt(float(vec2 * vec2))
    angl = acos(float(multi)/float(length))
    return angl

def calc_distance(v1,v2):
    dis = sqrt((v1-v2)*(v1-v2))
    return dis

residue_number_list=[45,48,49,52,53,56,57,60,61]

pyrosetta.init()
protein=pyrosetta.rosetta.core.import_pose.pose_from_file("orient_PP2_input.pdb")
scaffold=pyrosetta.rosetta.core.import_pose.pose_from_file("DHR7_test_protein_renumber.pdb")
rots = generate_canonical_rotamer_residues('TYR','helical')
print len(rots)
v=pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t()
v.append(('CE1','CE1'))
v.append(('CZ','CZ'))
v.append(('CE2','CE2'))
num = 0
t0=time.time()
final_ouput=[10.0]
'''
res1=[Vec(scaffold.residue(78).xyz("N")),Vec(scaffold.residue(78).xyz("CA")),Vec(scaffold.residue(78).xyz("C"))]
res2=[Vec(scaffold.residue(120).xyz("N")),Vec(scaffold.residue(120).xyz("CA")),Vec(scaffold.residue(120).xyz("C"))]
stub1 = stub(res1[0],res1[1],res1[2])
stub2 = stub(res2[0],res2[1],res2[2])
xform = stub2 * ~stub1
axis, ang, cen = xform.rotation_axis_center()
print axis.y
print ang
print cen.y
'''

for rot in rots:
    num = num+1
    print num
    print time.time() - t0
    rot.orient_onto_residue(protein.residue(4),v)
    #p=protein.clone()
    #p.replace_residue(4,rot,False)
    #p.dump_pdb('{}_test.pdb'.format(num))
    
    for i in residue_number_list:
        print i
        movable_coords = get_anchor_coordinates_from_pose(scaffold.residue(i).clone())
        fixed_coords = get_anchor_coordinates_from_pose(rot.clone())
        R,t = np_utils.rigid_transform_3D(movable_coords,fixed_coords)
        s=scaffold.clone()
        rotate_pose(s,numpy_to_rosetta(R))
        np_utils.translate_pose(s,numpy_to_rosetta(t.T))
        #s.dump_pdb("test_{}.pdb".format(i))
        res1=[Vec(s.residue(45).xyz("N")),Vec(s.residue(45).xyz("CA")),Vec(s.residue(45).xyz("C"))]
        res2=[Vec(s.residue(87).xyz("N")),Vec(s.residue(87).xyz("CA")),Vec(s.residue(87).xyz("C"))]
        stub1 = stub(res1[0],res1[1],res1[2])
        stub2 = stub(res2[0],res2[1],res2[2])
        xform = stub2 * ~stub1
        axis, ang, cen = xform.rotation_axis_center()
        angle_to_z = get_vector_angl(Vec(axis.x,axis.y,axis.z),Vec(0.0,0.0,1.0))
        dis_to_z = sqrt(cen.x*cen.x+cen.y*cen.y)
        if angle_to_z < final_ouput[0]:
            final_ouput=[angle_to_z,dis_to_z,i,rot.chi(1),rot.chi(2)]
            print final_ouput
        if angle_to_z < 0.1:
            s.dump_pdb("{}_{}_{}_{}_{}.pdb".format(angle_to_z,dis_to_z,i,rot.chi(1),rot.chi(2)))
print final_ouput







'''
v1=Vec(protein.residue(143).xyz('O')) - Vec(protein.residue(136).xyz('O'))
v2=Vec(protein.residue(102).xyz('O')) - Vec(protein.residue(95).xyz('O'))
v3=Vec(protein.residue(60).xyz('O')) - Vec(protein.residue(53).xyz('O'))
v4=Vec(protein.residue(15).xyz('O')) - Vec(protein.residue(8).xyz('O'))
print get_vector_angl(v1,Vec(0.0,0.0,1.0))
print get_vector_angl(v2,Vec(0.0,0.0,1.0))
print get_vector_angl(v3,Vec(0.0,0.0,1.0))
print get_vector_angl(v4,Vec(0.0,0.0,1.0))
print 'angle between CN and axis of helix'
v5=Vec(protein.residue(143).xyz('C')) - Vec(protein.residue(143).xyz('N'))
v6=Vec(protein.residue(142).xyz('C')) - Vec(protein.residue(142).xyz('N'))
v7=Vec(protein.residue(141).xyz('C')) - Vec(protein.residue(141).xyz('N'))
v8=Vec(protein.residue(140).xyz('C')) - Vec(protein.residue(140).xyz('N'))
v9=Vec(protein.residue(139).xyz('C')) - Vec(protein.residue(139).xyz('N'))
v10=Vec(protein.residue(138).xyz('C')) - Vec(protein.residue(138).xyz('N'))
v11=Vec(protein.residue(137).xyz('C')) - Vec(protein.residue(137).xyz('N'))
print get_vector_angl(v1,v5)
print get_vector_angl(v1,v6)
print get_vector_angl(v1,v7)
print get_vector_angl(v1,v8)
print get_vector_angl(v1,v9)
print get_vector_angl(v1,v10)
print get_vector_angl(v1,v11)
protein.dump_pdb("DHR7_test_protein_renumber.pdb")
'''

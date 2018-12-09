#!/home/sheffler/venv/david/bin/ipython
import pyrosetta
from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
from xyzMath import *
import math
import pyrosetta.toolbox.numpy_utils as np_utils
from numpy_utils import *
import random

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
    tryrot = TryRotamers(3, sf, 3, 0, True, solo_res=False, include_current=False )
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

def get_vector_angl(vec1,vec2):
    multi = vec1 * vec2
    length = sqrt(float(vec1 * vec1)) * sqrt(float(vec2 * vec2))
    angl = acos(float(multi)/float(length))
    return angl

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

def get_anchor_coordinates_from_pose(residue):
    bb_atoms = ['N','CA','C']
    coords = []
    for atom in bb_atoms:
        coords.append([residue.xyz(atom).x,residue.xyz(atom).y,residue.xyz(atom).z])
    return np.mat(coords)

def calc_distance(v1,v2):
    dis = sqrt((v1-v2)*(v1-v2))
    return dis

pyrosetta.init(" -extra_res_fa CAL.params CO3.params ")
p=pyrosetta.rosetta.core.import_pose.pose_from_file("final.pdb")
scaffold=pyrosetta.rosetta.core.import_pose.pose_from_file("mute_6k+1_6k+5.pdb")
surface=pyrosetta.rosetta.core.import_pose.pose_from_file("orient_final_surface.pdb")

'''
v1 = Vec(scaffold.residue(13).xyz('CA'))
v2 = Vec(scaffold.residue(13).xyz('HA'))
vec = v2 - v1

print vec
angle = get_vector_angl(vec, Vec(0.0,0.0,1.0))
print angle
'''
Vec_surface_1 = Vec(p.residue(81).xyz("Ca2p"))-Vec(p.residue(85).xyz("Ca2p"))
Vec_surface_2 = Vec(p.residue(62).xyz("Ca2p"))-Vec(p.residue(85).xyz("Ca2p"))
Vec_surface_normal = get_normal_vector(p.residue(81).xyz("Ca2p"),p.residue(85).xyz("Ca2p"),p.residue(62).xyz("Ca2p"))
#print get_vector_angl(vec,Vec_surface_normal)

v=pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t()
v.append(('OE1','OE1'))
v.append(('CD','CD'))
v.append(('OE2','OE2'))
rots = generate_canonical_rotamer_residues('GLU',(-118.6,118.3))
num = 0
angle_bias = [0.0,0.0,999.0]
for rot in rots:
    bias = []
    num = num + 1
    print num
    rot.orient_onto_residue(p.residue(139),v)
    p.replace_residue(139,rot,False)
    NC_vector = Vec(p.residue(139).xyz('C'))-Vec(p.residue(139).xyz('N'))
    CC_vector = Vec(p.residue(139).xyz('CA'))-Vec(p.residue(139).xyz('CB'))
    CH_vector = Vec(p.residue(139).xyz('CA'))-Vec(p.residue(139).xyz('HA'))
    bias.append(abs(get_vector_angl(NC_vector,Vec_surface_1)-1.528))
    bias.append(abs(get_vector_angl(CC_vector,Vec_surface_1) - 1.839))
    if (bias[0]+bias[1]) < angle_bias[2] and get_vector_angl(CH_vector,Vec_surface_1) < 1.57:
        angle_bias = [bias[0],bias[1],float(bias[0]+bias[1]),Vec(p.residue(139).xyz('C')),Vec(p.residue(139).xyz('N')),Vec(p.residue(139).xyz('CA')),p.chi(1,139),p.chi(2,139),p.chi(3,139),get_vector_angl(NC_vector,Vec_surface_1),get_vector_angl(CC_vector,Vec_surface_normal)]
    #print get_vector_angl(NC_vector,Vec_surface_1)
    #print get_vector_angl(CH_vector,Vec_surface_normal)
    #print p.chi(1,139)
    #print p.chi(2,139)
    #print p.chi(3,139)
    #p.dump_pdb('{}.pdb'.format(num))
print angle_bias


p.dump_pdb("fifi.pdb")
q=p.residue(139).clone()
q.set_chi(1,angle_bias[6])
q.set_chi(2,angle_bias[7])
q.set_chi(3,angle_bias[8])
q.orient_onto_residue(p.residue(139),v)
p.replace_residue(139,q,False)
p.dump_pdb("chosedone.pdb")
mute = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
mute.set_res_name('GLU')
res_nums = [x for x in range(2,scaffold.size()) if x%18 == 13]
for i in res_nums:
    mute.set_target(i)
    mute.apply(scaffold)
    scaffold.set_chi(1,i,angle_bias[6])
    scaffold.set_chi(2,i,angle_bias[7])
    scaffold.set_chi(3,i,angle_bias[8])
movable_coords = get_anchor_coordinates_from_pose(scaffold.residue(13).clone())
fixed_coords = get_anchor_coordinates_from_pose(p.residue(139).clone())
R,t = np_utils.rigid_transform_3D(movable_coords,fixed_coords)
rotate_pose(scaffold,numpy_to_rosetta(R))
np_utils.translate_pose(scaffold,numpy_to_rosetta(t.T))
p.delete_polymer_residue(139)
p.append_pose_by_jump(scaffold,138)
p.dump_pdb("final_scaffold_onsurface_3.pdb")
'''
virtual_OE2 = Vec(p.residue(139).xyz('OE2')) + 9.8*Vec_surface_1/sqrt(Vec_surface_1*Vec_surface_1)
virtual_CD = Vec(p.residue(139).xyz('CD')) + 9.8*Vec_surface_1/sqrt(Vec_surface_1*Vec_surface_1)
virtual_OE1 = Vec(p.residue(139).xyz('OE1')) + 9.8*Vec_surface_1/sqrt(Vec_surface_1*Vec_surface_1)
mute = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
mute.set_res_name('GLU')
res_nums = [x for x in range(2,scaffold.size()) if x%18 == 13]
for i in res_nums:
    mute.set_target(i)
    mute.apply(scaffold)
    scaffold.set_chi(1,i,angle_bias[6])
    scaffold.set_chi(2,i,angle_bias[7])
    scaffold.set_chi(3,i,angle_bias[8])
scaffold.dump_pdb("set_GLU.pdb")
random_center_1=angle_bias[6]
random_center_2=angle_bias[7]
random_center_3=angle_bias[8]
movable_coords = get_anchor_coordinates_from_pose(scaffold.residue(13).clone())
fixed_coords = get_anchor_coordinates_from_pose(p.residue(139).clone())
R,t = np_utils.rigid_transform_3D(movable_coords,fixed_coords)
rotate_pose(scaffold,numpy_to_rosetta(R))
np_utils.translate_pose(scaffold,numpy_to_rosetta(t.T))
rmsd=calc_distance(Vec(scaffold.residue(49).xyz('OE2')),virtual_OE2)+calc_distance(Vec(scaffold.residue(49).xyz('CD')),virtual_CD)+calc_distance(Vec(scaffold.residue(49).xyz('OE1')),virtual_OE1)
print rmsd
for i in range(1000):
    chi_1 = random.uniform(-180.0,180.0)
    chi_2 = random.uniform(random_center_2-20.0,random_center_2+20.0)
    #chi_3 = random.uniform(random_center_3-10.0,random_center_3+10.0)
    for j in res_nums:
        scaffold.set_chi(1,j,chi_1)
        scaffold.set_chi(2,j,chi_2)
        #scaffold.set_chi(3,j,chi_3)
    movable_coords = get_anchor_coordinates_from_pose(scaffold.residue(13).clone())
    fixed_coords = get_anchor_coordinates_from_pose(p.residue(139).clone())
    R,t = np_utils.rigid_transform_3D(movable_coords,fixed_coords)
    rotate_pose(scaffold,numpy_to_rosetta(R))
    np_utils.translate_pose(scaffold,numpy_to_rosetta(t.T))
    rmsd_temp=calc_distance(Vec(scaffold.residue(49).xyz('OE2')),virtual_OE2)+calc_distance(Vec(scaffold.residue(49).xyz('CD')),virtual_CD)+calc_distance(Vec(scaffold.residue(49).xyz('OE1')),virtual_OE1)
    #print rmsd_temp
    if rmsd_temp < rmsd:
        random_center_1=chi_1
        random_center_2=chi_2
        #random_center_3=chi_3
        print rmsd_temp
        rmsd = rmsd_temp

p.delete_polymer_residue(139)
p.append_pose_by_jump(scaffold,138)
p.dump_pdb("final_scaffold_onsurface.pdb")
'''


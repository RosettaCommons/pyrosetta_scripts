#!/home/sheffler/venv/david/bin/python
import pyrosetta
import argparse
from math import *
from xyzMath import *
import random
import collections
from rosetta import *
from rosetta.numeric import xyzVector_double_t as V3
from rosetta.numeric import xyzMatrix_double_t as M3
import pickle
from rosetta.protocols.protein_interface_design.movers import TryRotamers
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from hash_subclass import *
import time

pro_C_H_pair = [(7,14)]
pdb_name_list = ["origin_complex.pdb"]

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

def rotate_around_axis_for_tyrosine(angle,p):
    v_normal = get_normal_vector(p.xyz(pyrosetta.rosetta.core.id.AtomID(7,4)),p.xyz(pyrosetta.rosetta.core.id.AtomID(11,4)),p.xyz(pyrosetta.rosetta.core.id.AtomID(10,4)))
    arom_center = p.xyz(pyrosetta.rosetta.core.id.AtomID(7,4)) + p.xyz(pyrosetta.rosetta.core.id.AtomID(11,4)) +p.xyz(pyrosetta.rosetta.core.id.AtomID(10,4))
    arom_center_vec = Vec(float(arom_center[0])/3.0, float(arom_center[1])/3.0, float(arom_center[2])/3.0)
    virtual_atom = arom_center_vec + v_normal
    prot_bb=V3(virtual_atom.to_rosetta()),V3(arom_center_vec.to_rosetta()),V3(p.residue(4).xyz('CZ'))
    frame = pyrosetta.rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
    virtual_local = Vec(frame.global2local(virtual_atom.to_rosetta()))
    Center_local = Vec(frame.global2local(arom_center_vec.to_rosetta()))
    normal_axis = virtual_local - Center_local
    t= rotation_matrix_degrees(normal_axis,angle)
    for ia in range(1,p.clone().residue(4).natoms()+1):
        aid = pyrosetta.rosetta.core.id.AtomID(ia,4)
        atom_local = Vec(frame.global2local(p.xyz(aid)))
        atom_local_new = t*atom_local
        atom_global_new = frame.local2global(atom_local_new.to_rosetta())
        p.set_xyz(aid,atom_global_new)
    return p

def get_vector_angl(vec1,vec2):
    multi = abs(vec1 * vec2)
    length = sqrt(float(vec1 * vec1)) * sqrt(float(vec2 * vec2))
    angl = acos(float(multi)/float(length))
    return angl

def get_ala_asn_from_pdb(p):
    q=move_jump(p.clone())
    return q
#set up the move_jump
def move_jump(p):
    chains=p.split_by_chain()
    Surf=chains[1]
    Resi=chains[2]
    #change here
    Surf.append_pose_by_jump(Resi,2,'N','CG')
    return Surf
#set up score function
def sc_energy_fxn2():
    sf = pyrosetta.get_score_function()
    #print sf
    sf.set_weight(core.scoring.fa_sol, 0.1)
    sf.set_weight(core.scoring.fa_rep, 0.55)
    #sf.set_weight(core.scoring.fa_elec, 1.0)
    #sf.set_weight(core.scoring.hbond_sc, 1.0)
    sf.set_weight(core.scoring.hbond_bb_sc, 1.0)
    #sf.set_weight(core.scoring.atom_pair_constraint, 1.0)
    #sf.set_weight(core.scoring.angle_constraint, 1.0)
    #sf.set_weight(core.scoring.dihedral_constraint, 1.0)
    return sf

def get_constraint_parameters(p,pairs):
    #pro_C_H_pair = [(2,9),(5,10),(5,11),(6,12),(6,13),(7,14),(7,15)]
    #for pairs in pro_C_H_pair:
    C_ID,H_ID = pairs
    v_normal = get_normal_vector(p.residue(4).xyz('CG'),p.residue(4).xyz('CE1'),p.residue(4).xyz('CE2'))
    v_H = Vec(p.xyz(pyrosetta.rosetta.core.id.AtomID(H_ID,2)))
    v_C = Vec(p.xyz(pyrosetta.rosetta.core.id.AtomID(C_ID,2)))
    v_CH = v_C - v_H
    ang = get_vector_angl(v_normal , v_CH)
    arom_center = p.residue(4).xyz('CG') + p.residue(4).xyz('CE1') + p.residue(4).xyz('CE2')
    #arom_center = p.xyz(pyrosetta.rosetta.core.id.AtomID(7,4)) + p.xyz(pyrosetta.rosetta.core.id.AtomID(11,4)) +p.xyz(pyrosetta.rosetta.core.id.AtomID(10,4))
    arom_center_vec = Vec(float(arom_center[0])/3.0, float(arom_center[1])/3.0, float(arom_center[2])/3.0)
    Hcenter_vec = v_H - arom_center_vec
    Ccenter_vec = v_C - arom_center_vec
    Hcenter_distance = sqrt(float(Hcenter_vec * Hcenter_vec))
    Ccenter_distance = sqrt(float(Ccenter_vec * Ccenter_vec))
    H_ang = get_vector_angl(v_normal , Hcenter_vec)
    H_projection_distance = float(Hcenter_distance * sin(H_ang))
    return (Hcenter_distance,H_projection_distance,ang,Ccenter_distance)

def get_energy_score_PP2(Hcenter_distance,H_projection_distance,ang,Ccenter_distance):
    if Ccenter_distance <= Hcenter_distance:
        extra_score = 999.0
    else:
        if float(Hcenter_distance)>4.0:
            Hcenter_score = 50*(float(Hcenter_distance)-4.0)*(float(Hcenter_distance)-4.0)+1.62
        else:
            Hcenter_score = 2*(float(Hcenter_distance)-3.1)*(float(Hcenter_distance)-3.1)
        #print float(Hcenter_distance)
        if H_projection_distance > 2.0:
            H_projection_score = 5*(H_projection_distance-2.0)+1.85
        else:
            H_projection_score = abs(0.15 - H_projection_distance)
        #print H_projection_distance
        if float(ang) > 1.4:
            ang_score = 50*(float(ang)-1.4)*(float(ang)-1.4)+0.98
        else :
            ang_score = 2*(float(ang)-0.7)*(float(ang)-0.7)
        #print ang
        extra_score = Hcenter_score + H_projection_score + ang_score
    return extra_score

def rotation_around_axis(angle, p):
    prot_bb=V3(p.residue(4).xyz('CZ')),V3(p.residue(4).xyz('OH')),V3(p.residue(4).xyz('HH'))
    frame = pyrosetta.rosetta.core.kinematics.Stub(prot_bb[0],prot_bb[1],prot_bb[2])
    Hatom_local = Vec(frame.global2local(p.residue(4).xyz('HH')))
    #print Hatom_local
    Catom_local = Vec(frame.global2local(p.residue(4).xyz('CZ')))
    #print Catom_local
    Oatom_local = Vec(frame.global2local(p.residue(4).xyz('OH')))
    #print Oatom_local
    v1 = Oatom_local
    v2 = Catom_local
    CO_vec = v1 - v2
    #print CO_vec
    t=rotation_matrix_degrees(CO_vec,angle)
    Hatom_local_new = t*Hatom_local
    Hatom_global_new = frame.local2global(Hatom_local_new.to_rosetta())
    aid = pyrosetta.rosetta.core.id.AtomID(24,4)
    p.set_xyz(aid, Hatom_global_new)
    return p


#set up the parameters
options={}
#options['input_pdb']='PP2_CAHA.pdb'
options['energy_function']='-empty'
options['constraint_file']=''
options['extra_res_fa']=''
options['jump']=1
options['Gaus_trans_mag']=float(0.5)
options['Gaus_rot_mag']=float(3)
options['Energy_cut']=6
options['Resi_num']=4
options['Resi_name']='TYR'
options['hashtable_name']='PP2_hashtable_1'
options['target_phi_psi_list'] = [(-90.0,-2.0)]
options['orient_atoms'] = ['CE1','CZ','CE2']
options['create_ray_atoms_num']=[2,2,2,2]
options['create_ray_atoms_name']=['C','O','CA','N']
options['chi_angle_num']=2


if __name__ == '__main__':
    #pairs = (2,9)
    t0=time.time()
    print '[ATTENTION] Please make sure you have set up move_jump, all the parameters and the score function in this script!'
    pdb_jump = int(options['jump'])
    opts=[]
    opts.append(options['energy_function'])
    if options['extra_res_fa']:
        opts.append('-extra_res_fa')
        opts.append(options['extra_res_fa'])
    if options['constraint_file']:
        opts.append('-cst_file')
        opts.append(options['constraint_file'])
    print options
    print opts
    pyrosetta.init(extra_options=' '.join(opts))
    resl = 1.0
    lever = 3.0
    hc = make_hash_store_rot(resl , lever, options['target_phi_psi_list'],options['orient_atoms'])
    for m in range(0,1):
        print 'number {} structure begins!'.format(m)
        pairs = pro_C_H_pair[m]
        options['input_pdb']=pdb_name_list[m]
        for phi in range(-80,-50,5):
            for psi in range(145,175,5):
                print time.time()-t0
                p=pyrosetta.rosetta.core.import_pose.pose_from_file(options['input_pdb'])
                gn=get_ala_asn_from_pdb(p)
                gn.set_phi(2,float(phi))
                gn.set_psi(2,float(psi))
                gn.dump_pdb('origin_complex.pdb')
                if options['constraint_file']:
                    pyrosetta.rosetta.core.scoring.constraints.add_constraints_from_cmdline_to_pose(gn)
                j=gn.jump(pdb_jump)
                #sample_resi = gn.residue(int(options['Resi_num']))
                print '[ATTENTION] Please check residue numbers in the constraint file!'
                sf_sys = sc_energy_fxn2()
                old_trans = V3(j.get_translation())
                old_rot = M3(j.get_rotation())
                low_rot=old_rot
                low_trans=old_trans
                parameters=get_constraint_parameters(gn,pairs)
                old_E=sf_sys(gn) + float(get_energy_score_PP2(*parameters))
                #old_E=sf_sys(gn)
                low_E=old_E
                print(old_E)
                #print(sf_sys)
                T = 0.6
                total_cycle = int(10000)
                n_accept = 0
                for i in range(total_cycle):
                    rand=random.random()
                    k = 1
                    if rand <0.5 :
                        k = -1
                    j.gaussian_move(k,options['Gaus_trans_mag'],options['Gaus_rot_mag'])
                    gn.set_jump(pdb_jump,j)
                    parameters=get_constraint_parameters(gn,pairs)
                    new_E = sf_sys(gn) + float(get_energy_score_PP2(*parameters))
                    gn_copy = gn
                    #sample_resi = gn.residue(int(options['Resi_num']))
                    #name='test_gaussian_' + str(i) + '.pdb'
                    #gn.dump_pdb(name)
                    for rot in range(0,360,10):
                        rot_ang = float(rot)
                        gn_rot = rotation_around_axis(rot_ang, gn_copy.clone())
                        #name = "rotate_{}.pdb".format(rot)
                        #gn_rot.dump_pdb(name)
                        #print 'nono'
                        if sf_sys(gn_rot) + float(get_energy_score_PP2(*parameters))< new_E:
                            #print 'yesyes'
                            new_E = sf_sys(gn_rot) + float(get_energy_score_PP2(*parameters))
                            gn = gn_rot
                            #print new_E

                    delta = new_E - old_E
                    accept = 1
                    if delta > 0 :
                        if math.exp(-delta/T) <random.random():
                            accept = 0
                    if accept :
                        #name='accept_gaussian_' + str(i) + '.pdb'
                        #gn.dump_pdb(name)
                        print new_E
                        old_E = new_E
                        old_trans = V3(j.get_translation())
                        old_rot = M3(j.get_rotation())
                        if low_E > new_E:
                            low_E = new_E
                            low_rot = old_rot
                            low_trans = old_trans

                        if new_E < float(options['Energy_cut']) :
                            #put it into hash
                            bb_rays = []
                            for num in range(4):
                                bb_rays.append(V3(gn.residue(int(options['create_ray_atoms_num'][num])).xyz(options['create_ray_atoms_name'][num])))
                            #sample_resi = gn.residue(int(options['Resi_num']))
                            Mark = hc.update_table_chi(bb_rays,gn.residue(int(options['Resi_num'])),options['Resi_name'],options['chi_angle_num'])
                            print str(Mark)
                            n_accept = n_accept + 1
                            if str(Mark) == 'True':
                                break;
                    else:
                        j.set_rotation(old_rot)
                        j.set_translation(old_trans)
                        gn.set_jump(pdb_jump,j)

        print('low energy: %s'%low_E)
        j.set_rotation(low_rot)
        j.set_translation(low_trans)
        gn.set_jump(pdb_jump,j)
        rot1 = make_hash(resl,lever,options['target_phi_psi_list'],options['orient_atoms']).orient_rots(gn.residue(int(options['Resi_num'])), options['Resi_name'])[1]
        gn.dump_pdb('finalpre_{}.pdb'.format(m))
        gn.replace_residue(int(options['Resi_num']),rot1,False)
        gn.dump_pdb('final_{}.pdb'.format(m))
    print "Now build hashtable!"
    pickle.dump(hc.dd, open(options['hashtable_name'], "wb"))
    print time.time()-t0



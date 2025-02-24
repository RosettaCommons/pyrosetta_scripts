from pdb_utils_noclass import *
from pyrosetta.toolbox import pose_from_rcsb
import random
from math import *
from rif.legacy.xyzMath import *
from rosetta.protocols.sic_dock import Rose
from rosetta.core.id import AtomID
from rosetta.core.pose import xyzStripeHashPose, PoseCoordPickMode_BB, PoseCoordPickMode_CB
from rosetta.numeric import xyzVector_float_t
from rif import stripe_index_3d, rcl
from rosetta.numeric import xyzVector_float_t
from rosetta.core.kinematics import Stub


class PoseResnumIndex(object):
    def __init__(self, pose, clash_dis=3.2, nbr_dis=8.5, nbr_atom='CA'):
        bb_atoms = rcl.atoms(pose, 'bbheavy')
        self.clash_check = stripe_index_3d(clash_dis, bb_atoms)
        cens, resnums = list(), list()
        for res in pose:
            if res.is_protein():
                cens.append(res.xyz('CA'))
                resnums.append(res.seqpos())
        self.resi_map = stripe_index_3d(nbr_dis, cens, resnums)
        self.nbr_atom = nbr_atom

    def clashes(self, pose):
        for res in pose:
            for ia in range(1, res.nheavyatoms() + 1):
                xyz = res.xyz(ia)
                if self.clash_check.neighbor_exists(xyz):
                    return True
        return False

    def contacting_resnums(self, pose):
        if self.clashes(pose):
            return None  # could alternately return empty set
        contacts = set()
        for res in pose:
            if res.has(self.nbr_atom):
                xyz = res.xyz(self.nbr_atom)
                nbrs = self.resi_map.neighbors(xyz)
                contacts.update(nbrs)
        return contacts


class ClashAndContact(object):
    def __init__(self, pose):
        self.clash_check = xyzStripeHashPose(pose, PoseCoordPickMode_BB, 3.2)
        self.cb_check = xyzStripeHashPose(pose, PoseCoordPickMode_CB, 8.0)

    def contact_score(self, pose):
        for ir in range(1, pose.size() + 1):
            for ia in range(1, pose.residue(ir).nheavyatoms() + 1):
                xyz = pose.xyz(AtomID(ia, ir))
                if self.clash_check.clash(xyzVector_float_t(xyz.x, xyz.y, xyz.z)):
                    return -1
        contacts = 0
        for ir in range(1, pose.size() + 1):
            xyz = pose.xyz(AtomID(4, ir))
            contacts += self.cb_check.nbcount(xyzVector_float_t(xyz.x, xyz.y, xyz.z))
        return contacts

def eval_contacts_pair(rep_pose,dock_gen,contact_cutoff):
   good_matches=[]
   contact_hist={}
   checker=PoseResnumIndex(rep_pose)
   ntries=0
   for pept_pose, deg, dist in dock_gen:
         ntries += 1
         contact_set = checker.contacting_resnums(pept_pose)

         if contact_set is None:
             contacts = 0
         else:
            contacts=len(contact_set)

         if contacts in list(contact_hist.keys()):
             contact_hist[contacts]=contact_hist[contacts]+1
         else:
             contact_hist[contacts]=1
#         print contacts
#         print complex[1][0:5], contacts
         if contacts > contact_cutoff:
             good_matches.append( (pept_pose,deg, dist  ,contacts, contact_set) )
   return good_matches, contact_hist, ntries

def eval_contacts(complexes):
   first=complexes[0][1][0:5]  #DHR name
   rep_pose=complexes[0][0][2]
   checker=ClashAndContact(rep_pose)

   for complex in complexes:
         name=complex[1][0:5]
         if name != first:
             rep_pose=complex[0][2]
             first=name
             checker=CLashAndContact(rep_pose)
         pept_pose=complex[0][1]

         contacts = checker.contact_score(pept_pose)
         print(complex[1][0:5], contacts)

def ran_range(frange):
    return random.uniform(-frange,+frange)


def rotation_axis_center(X):
        axis, ang = X.R.rotation_axis()

        # these points lie on a circle who's center is the center of rotation
        p0 = Vec(0, 0, 0)
        p1 = Vec(X * p0)
        p2 = Vec(X * p1)
        p1 -= axis * (p1 - p0).dot(axis)
        p2 -= axis * (p2 - p1).dot(axis)

        assert(abs((p1 - p0).dot(axis)) < 0.000001)
        assert(abs((p2 - p1).dot(axis)) < 0.000001)

        d = p1.length()

        if(d < 0.000001):
            return axis, ang, Vec(0, 0, 0)

        if abs(2.0 * math.tan(ang / 2.0)) < 0.001:
            return axis, ang, None

        l = d / (2.0 * math.tan(ang / 2.0))

        tocen = p1.normalized().cross(axis) * l
        assert(abs(tocen.length() - l) < 0.0001)

        # correct direction based on curvature
        if(tocen.dot(p2 - p1) < 0.0):
            tocen = -tocen

        cen = (p0 + p1) / 2.0 + tocen
        return axis, ang, cen


def get_helix_params(res1,res2):
    #print "DEBUG, res1:", res1
    #print "DEBUG, res2:", res2
    stub1 = stub(res1[0],res1[1],res1[2])
    stub2 = stub(res2[0],res2[1],res2[2])
    xform = stub2 * ~stub1
    #Some error here:
    #File "/home/sheffler/venv/david/local/lib/python2.7/site-packages/xyzMath.py", line 892, in rotation_axis_center
    #assert(abs((p1 - p0).dot(axis)) < 0.000001)
    #print "DEBUG: stub1:", stub1
    #print "DEBUG: stub2:", stub2
    #print "DEBUG: xform:", xform

    axis, ang, cen = xform.rotation_axis_center()

    translation_along_axis_of_rotation = axis.dot(xform.t)
    radius_from_CA_1 = projperp(axis, cen - res1[1] ).length()
    if translation_along_axis_of_rotation < 0.0:
        ang = ang * -1.
        translation_along_axis_of_rotation=translation_along_axis_of_rotation * -1.

    return translation_along_axis_of_rotation,radius_from_CA_1,ang




def get_complexes_from_list(c_file):
 p=Pose()
 complex_file= list(map(str.split,open(c_file,'r').readlines()))
 complexes=[]
 for line in complex_file:
     pdb=line[0]
     p=rosetta.core.import_pose.pose_from_file(pdb)
     chains=p.split_by_chain()
     complexes.append( (chains,pdb) )
 return complexes


def eval_clashes_and_contacts(complexes,ala_convert):
   first=complexes[0][1]
   rep_pose=complexes[0][0][2]
   if ala_convert:
     rep_pose.dump_pdb('repeat.pdb')
     rep_pose=convert_to_ala(rep_pose.clone())
     rep_pose.dump_pdb('repeat_ala.pdb')
   repeat_rose=Rose(rep_pose) # just take repeat protein from 1st complex
#   clashes=[]
#   contacts=[]
   for complex in complexes:
         name=complex[1]
         if name != first:
             rep_pose=complex[0][2]
             first=name
             if ala_convert:
                 rep_pose=convert_to_ala(rep_pose.clone())
         repeat_rose=Rose(rep_pose)
         pept_pose=complex[0][1]
         if ala_convert:
             pept_pose=convert_to_ala(pept_pose.clone())
         pept_rose=Rose(pept_pose)
      #   clashes.append(repeat_rose.clashes(pept_rose))
         clashes = repeat_rose.clashes(pept_rose)
       #  contacts.append(repeat_rose.contacts(pept_rose))
         contacts = repeat_rose.contacts(pept_rose)
         print((name,clashes,contacts))

#Deprecated, this function is too slow because it generates a efxn every time
def eval_rama(aa,phi,psi):
    sf = rosetta.core.scoring.Ramachandran()
    res_type = rosetta.core.chemical.aa_from_name(aa)
    score=sf.eval_rama_score_residue(res_type,phi,psi)
    return score

# have to wait for python 3 for this
#def calc_repeat_protein_params_ws(input_file, struct_dir, *, offset=0):
def calc_repeat_protein_params_ws(input_file, struct_dir, offset):
 verbose=1
 p=Pose()
 rep_file=list(map(str.split,open(input_file,'r').readlines()))
 helical_params=[]
 helical_params2=[]
 helical_params3=[]
 names=[]
 lengths=[]
 pdbs={}

 for line in rep_file:

     # struct_dir = "/work/baker/repeat_peptide/designs/"
     DHR_file=struct_dir+line[0]
     DHR=str.split(line[0],'.')[0]

     length=int(line[1])
     lengths.append(length)

     names.append(DHR)

     if verbose: print(DHR_file,DHR,length)
     ## p=rosetta.core.import_pose.pose_from_file(DHR_file)
     p=convert_to_ala(rosetta.core.import_pose.pose_from_file(DHR_file))

     #for helix param calculation, go from middle of 2nd repeat to middle of 3rd repeat
     half_repeat=int(length/2)
     repeat1_start = length+half_repeat + offset
     repeat1_stop = repeat1_start + length-1
     repeat2_start = 2*length+half_repeat + offset
     repeat2_stop = repeat2_start+length-1

     repeat1_sel = rosetta.utility.vector1_unsigned_long()
     repeat2_sel = rosetta.utility.vector1_unsigned_long()
     for i in range(repeat1_start, repeat1_stop + 1):
        repeat1_sel.append(i)
     for i in range(repeat2_start, repeat2_stop + 1):
        repeat2_sel.append(i)

     ft = rosetta.core.kinematics.FoldTree(repeat1_stop - repeat1_start + 1)
     repeat1 = rosetta.core.pose.Pose()
     rosetta.core.pose.create_subpose(p, repeat1_sel, ft, repeat1)
     repeat2 = rosetta.core.pose.Pose()
     rosetta.core.pose.create_subpose(p, repeat2_sel, ft, repeat2)
     if verbose: print(p.size(),repeat1.size(),repeat2.size())
     # repeat1 and repeat2 are poses with adjacent repeat
    # segments, with identical length
    # now create an exact copy of repeat1, then superimpose it onto repeat2
    # this is basically a way to get the superposition xform, which is not
    # readily available from rosetta (that I am aware of)
     repeat1_onto_2 = rosetta.core.pose.Pose(repeat1)
     rms = rosetta.core.scoring.calpha_superimpose_pose(repeat1_onto_2, repeat2)
     if verbose: print('rms is', rms)
     res1 = [Vec(repeat1.residue(1).xyz('N')),  Vec(repeat1.residue(1).xyz('CA')), Vec(repeat1.residue(1).xyz('C'))]
     res2 = [Vec(repeat1_onto_2.residue(1).xyz('N')),Vec(repeat1_onto_2.residue(1).xyz('CA')), Vec(repeat1_onto_2.residue(1).xyz('C'))]
     trans, radius, ang = get_helix_params(res1,res2)
     helical_params.append( (trans, radius, ang) )
     p=center_on_z_axis(res1,res2,p.clone())

    #  p.dump_pdb('%s_tf.pdb'%DHR)
     pdbs[DHR]=p.clone()

 print('name repeat_length trans radius angle arc_length')
 helical_arcs={}
 for i in range(len(helical_params)):
    arc_length=sqrt( (helical_params[i][0])**2 +( ((helical_params[i][1]) * sin(helical_params[i][2] ))**2 ) )
    print(('%s %s   %.2f %.2f %.2f  %.2f'%(names[i],lengths[i],helical_params[i][0],helical_params[i][1],helical_params[i][2],arc_length)))
    helical_arcs[names[i]]='%s %s   %.2f %.2f %.2f  %.2f'%(names[i],lengths[i],helical_params[i][0],helical_params[i][1],helical_params[i][2],arc_length)
 return helical_params, helical_arcs, names, lengths, pdbs



def calc_repeat_protein_params():
 p=Pose()
 rep_file=list(map(str.split,open('repeat.list','r').readlines()))
 helical_params=[]
 helical_params2=[]
 helical_params3=[]
 names=[]
 lengths=[]
 first=24
 for line in rep_file:
     DHR=line[0]
     length=int(line[1])
     lengths.append(length)
     names.append(DHR)
     p=rosetta.core.import_pose.pose_from_file('../designs/%s'%DHR)
     res1=[Vec(p.residue(first).xyz("N")),Vec(p.residue(first).xyz("CA")),Vec(p.residue(first).xyz("C"))]
     res2=[Vec(p.residue(length+first).xyz("N")),Vec(p.residue(length+first).xyz("CA")),Vec(p.residue(length+first).xyz("C"))]
     res3=[Vec(p.residue(first+2*length).xyz("N")),Vec(p.residue(first+2*length).xyz("CA")),Vec(p.residue(first+2*length).xyz("C"))]
     res4=[Vec(p.residue(first+3*length).xyz("N")),Vec(p.residue(first+3*length).xyz("CA")),Vec(p.residue(first+3*length).xyz("C"))]
     trans1, radius1, ang1 = get_helix_params(res1,res2)
     trans2, radius2, ang2 = get_helix_params(res2,res3)
     trans3, radius3, ang3 = get_helix_params(res3,res4)
     helical_params.append( (trans1, radius1, ang1) )
     helical_params2.append( (trans2, radius2, ang2) )
     helical_params3.append( (trans3, radius3, ang3) )
     p=center_on_z_axis(res2,res3,p)

     p.dump_pdb('%s_tf.pdb'%DHR[0:base])

 return helical_params, helical_params2, helical_params3,names,lengths

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
            aid = rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(pose.xyz(aid))
            newxyz = transform * oldxyz
            pose.set_xyz(aid, newxyz.to_rosetta())

    return center_com_z(pose.clone())

def rotate_around_z(axis,angle,pose):
#    deg=2.*math.pi/360.
    r_axis=axis.to_rosetta()
#    transform=rosetta.numeric.rotation_matrix_double_t(r_axis,angle*deg)
#    v_axis=Vec(axis)
    t=rotation_matrix_degrees(Vec(0,0,1),angle)
    for ir in range(1, pose.size() + 1):
        for ia in range(1, pose.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(pose.xyz(aid))
            newxyz = t * oldxyz
            pose.set_xyz(aid, newxyz.to_rosetta())
    return pose


def translate_along_z(dist,pose):
    for ir in range(1, pose.size() + 1):
        for ia in range(1, pose.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia, ir)
            newxyz = Vec(pose.xyz(aid)) + Vec(0.,0.,dist)
            pose.set_xyz(aid, newxyz.to_rosetta())
    return pose

def convert_to_ala(input_pose):
    pose = input_pose.clone()
    chem_manager = rosetta.core.chemical.ChemicalManager
    rts = chem_manager.get_instance().residue_type_set("fa_standard")
    # print(rts)
    rfactory = rosetta.core.conformation.ResidueFactory
    res = rfactory.create_residue(rts.name_map('ALA'))
    for ir in range(1,pose.size() +1):
        pose.replace_residue(ir,res,True)
    return pose

def center_com_z(p):
   z_tot=0.
   n=0
   for ir in range(1,p.size()+1):
        for ia in range(1,p.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia, ir)
            z = p.xyz(aid).z
            n=n+1
            z_tot=z_tot+z
   z_vec=Vec(0.,0.,z_tot/float(n))
#   print('z offset %s'%z_vec)
   for ir in range(1, p.size() + 1):
        for ia in range(1, p.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia, ir)
            newxyz = Vec(p.xyz(aid))-z_vec
            p.set_xyz(aid, newxyz.to_rosetta())
   return p

def test():
    p=rosetta.core.import_pose.pose_from_file("DHR7.pdb")
    res1=[Vec(p.residue(50).xyz("N")),Vec(p.residue(50).xyz("CA")),Vec(p.residue(50).xyz("C"))]
    res2=[Vec(p.residue(92).xyz("N")),Vec(p.residue(92).xyz("CA")),Vec(p.residue(92).xyz("C"))]
    res3=[Vec(p.residue(134).xyz("N")),Vec(p.residue(134).xyz("CA")),Vec(p.residue(134).xyz("C"))]
    new_pose=center_on_z_axis(res2,res3,p.clone())
    p.dump_pdb('DHR7_tf.pdb')

def get_DHR7_and_arm_parameters():
    p=rosetta.core.import_pose.pose_from_file("DHR7.pdb")
    res1=[Vec(p.residue(50).xyz("N")),Vec(p.residue(50).xyz("CA")),Vec(p.residue(50).xyz("C"))]
    res2=[Vec(p.residue(92).xyz("N")),Vec(p.residue(92).xyz("CA")),Vec(p.residue(92).xyz("C"))]
    res3=[Vec(p.residue(134).xyz("N")),Vec(p.residue(134).xyz("CA")),Vec(p.residue(134).xyz("C"))]
    print('DHR7 params ')
    get_helix_params(res1,res2)
    get_helix_params(res2,res3)

    p=rosetta.core.import_pose.pose_from_file("ARM_pept.pdb")
    #nres=ARM_pept.size()
    res1=[Vec(p.residue(3).xyz("N")),Vec(p.residue(3).xyz("CA")),Vec(p.residue(3).xyz("C"))]
    res2=[Vec(p.residue(5).xyz("N")),Vec(p.residue(5).xyz("CA")),Vec(p.residue(5).xyz("C"))]
    print('ARM peptide params ')
    get_helix_params(res1,res2)

    p=rosetta.core.import_pose.pose_from_file("5AEI_A.pdb")
    print((p.residue(106).name(),p.residue(107).name(),p.residue(149).name(),p.residue(191).name()))
    res1=[Vec(p.residue(107).xyz("N")),Vec(p.residue(107).xyz("CA")),Vec(p.residue(107).xyz("C"))]
    res2=[Vec(p.residue(149).xyz("N")),Vec(p.residue(149).xyz("CA")),Vec(p.residue(149).xyz("C"))]
    res3=[Vec(p.residue(191).xyz("N")),Vec(p.residue(191).xyz("CA")),Vec(p.residue(191).xyz("C"))]
    print('ARM protein params ')
    get_helix_params(res1,res2)
    get_helix_params(res2,res3)
    return()

from rosetta.numeric import xyzVector_double_t as V3
from rosetta.numeric import xyzMatrix_double_t as M3
import rosetta.core.pack.rotamer_set
from rosetta.protocols.protein_interface_design.movers import TryRotamers

def generate_canonical_residue(residue_name3):
    work_pose = Pose()
    rosetta.core.pose.make_pose_from_sequence (work_pose, "X[%s]" % residue_name3, "fa_standard", auto_termini=False)
    work_residue = rosetta.core.conformation.Residue( work_pose.residue(1) )

    ca_loc = work_residue.xyz("CA")

    for a in range( work_residue.natoms() ):
        work_residue.set_xyz( a + 1, work_residue.xyz(a + 1) - ca_loc)

    return work_residue

def generate_canonical_rotamer_residues( residue_name3, target_phi_psi ):
    canonical_phi_psi = {
        "helical" : (-66.0, -42.0),
        "sheet"   : (-126.0, 124.0)
        }
    canonical_residue = generate_canonical_residue(residue_name3)
    test_sequence = "AAX[%s]AA" % residue_name3

    if target_phi_psi in canonical_phi_psi:
        target_phi, target_psi = canonical_phi_psi[ target_phi_psi ]
    else:
        target_phi, target_psi = target_phi_psi

    sf = pyrosetta.get_score_function()
    tryrot = TryRotamers(3, sf, 0, 0, True, solo_res=False, include_current=False ) #  3rd argument is rotamer explosion level.  2 means ex1 ex2

    test_pose = Pose()
    rosetta.core.pose.make_pose_from_sequence( test_pose, test_sequence, "fa_standard" )
    for i in range(1, test_pose.size() + 1):
        test_pose.set_psi(i, target_psi)
        test_pose.set_phi(i, target_phi)
        test_pose.set_omega(i, 180)

    tryrot.setup_rotamer_set( test_pose )

    rotamer_set = tryrot.rotamer_set()
    rotamers = [rotamer_set.rotamer(i).clone() for i in xrange(1, rotamer_set.num_rotamers() + 1)]

    for r in rotamers:
        r.orient_onto_residue( canonical_residue )
        r.seqpos(1)

    return rotamers


def get_richardson_rot_data():
 rr = pd.DataFrame(  # rr for richardson rotamers
    [
        ('asn', 'p10', 60, 62, -10, 8, 20, -90, 0),
        ('asn', 'p30', 100, 62, 30, 6, 20, 0, 90),
        ('asn', 't20', 100, -174, -20, 5, 20, -120, 0),
        ('asn', 't30', 228, -177, 30, 14, 15, 0, 80),
        ('asn', 'm20', 580, -65, -20, 10, 20, -60, 10),
        ('asn', 'm80', 118, -65, -75, 10, 10, -100, -60),
        ('asn', 'm120', 58, -65, 120, 10, 18, 60, 100),
    ], columns='aa rot count chi1 chi2 sd1 sd2 lb2 ub2'.split()
 )
 print(rr)
 assert rr.shape[0] == 7
 assert rr.shape[1] == 9
 #   print(rr.chi1[4])
 print("=== just listed rotamers ===")
 for i in range(rr.shape[0]):
        print('chi1/2', rr.chi1[i], rr.chi2[i])
 print("=== chi2 every 20 degrees ===")
 for i in range(rr.shape[0]):
        for chi2 in np.arange(rr.lb2[i], rr.ub2[i], 20.0):
            print('chi1/2', rr.chi1[i], chi2)
 return rr


def orient_rots(base,rots):
    #don't need new_rots and input rots are transformed in place
    new_rots=[]
    v=pyrosetta.rosetta.utility.vector1_std_pair_std_string_std_string_t()
    v.append( ('OD1', 'OD1') )
    v.append( ('ND2', 'ND2') )
    v.append( ('CG', 'CG') )
    for rot in rots:
        new_rots.append(rot.orient_onto_residue(base,v))
       
    return new_rots

def test_rot_generation(gn):
## test rot generation by superimposing on bb and dumping pdb
#x=rots[1].xyz('ND2')
 rots= generate_canonical_rotamer_residues("ASN","helical")
 new_rots=orient_rots(gn.residue(4),rots)  # test only.  orient rots on ASN from pdb
#x=rots[1].xyz('ND2')
 rot_pose=Pose()
 for rot in rots:
    print(rot.xyz('ND2'))
    rot_pose.append_residue_by_jump(rot,1,"","",True)
 rot_pose.dump_pdb('rots.pdb')
# end test

if name == "__main__":
#need to open pdb first
  test_rot_generation()

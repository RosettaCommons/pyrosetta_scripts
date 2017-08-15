#!/bin/env/python
"""
Based on PyRosetta's dna_interface.py

This script takes the structure files for PX DNA and a single head of
T7 endonuclease I and, with a given scissile phosphate will dock the two together
using PyRosetta's high-resolution docking.  This functions by first superimposing
PX DNA such that the scissile phosphate lines up with the scissile phosphate of
the Holliday Junction in the original T7 structure (2PFJ).  Next constraints are
generated to maintain this proper positioning followed by high resolution docking
which consists of small rigid-body perturbations, sidechain packing, and minimization.  This uses the "dna" score function.

Author: Benjamin D Walcott
    based on script by Evan H. Baugh
        based on an original script by Sid Chaudhury
        revised and motivated by Robert Schleif
    Updated by Boon Uranukul, 6/9/12
    Simplified special constant seed initialization ~ Labonte

References:
    J. J. Gray, "High-resolution protein-protein docking," Curr. Opinions in
        Struct. Bio. 16 (2) 183-193 (2006).

"""

#Imports
import optparse
#test for PyRosetta 4
PR4 = True
try:
    from pyrosetta import *
except ImportError:
    PR4 = False
from rosetta import *
init()
import rosetta.protocols.simple_moves
import os
import Bio
import Bio.PDB
import Bio.PDB.Structure
import Bio.PDB.Model
import shutil

def sample_dna_interface(pdb_filename, partners,
        jobs = 1, job_output = 'dna_output',cst_file = '',use_pymol = False):
    """
    This is more or less unchanged from dna_interface.py with the exception of
    the addition of constraints and the use_pymol parameter. - BDW

    Performs DNA-protein docking using Rosetta fullatom docking (DockingHighRes)
        on the DNA-protein complex in  <pdb_filename>  using the relative chain
        <partners>  .
        <jobs>  trajectories are performed with output structures named
        <job_output>_(job#).pdb.

    """
    # 1. creates a pose from the desired PDB file
    pose = Pose()
    #This function has changed in PyRosetta 4.  It still has pose_from_pdb for
    #backwards compatibility, but I'd prefer to use the new function.
    if PR4:
        pose_from_file(pose, pdb_filename)
    else:
        pose_from_pdb(pose,pdb_filename)

    # 2. setup the docking FoldTree
    # using this method, the jump number 1 is automatically set to be the
    #    inter-body jump
    dock_jump = 1
    # the exposed method setup_foldtree takes an input pose and sets its
    #    FoldTree to have jump 1 represent the relation between the two docking
    #    partners, the jump points are the residues closest to the centers of
    #    geometry for each partner with a cutpoint at the end of the chain,
    # the second argument is a string specifying the relative chain orientation
    #    such as "A_B" of "LH_A", ONLY TWO BODY DOCKING is supported and the
    #    partners MUST have different chain IDs and be in the same pose (the
    #    same PDB), additional chains can be grouped with one of the partners,
    #    the "_" character specifies which bodies are separated
    # the third argument...is currently unsupported but must be set (it is
    #    supposed to specify which jumps are movable, to support multibody
    #    docking...but Rosetta doesn't currently)
    # the FoldTrees setup by this method are for TWO BODY docking ONLY!
    protocols.docking.setup_foldtree(pose, partners, Vector1( [dock_jump]))

    # 3. create a copy of the pose for testing
    test_pose = Pose()
    test_pose.assign(pose)

    #Add constraints from a constraints file
    if cst_file != '':
        cstMover = protocols.simple_moves.ConstraintSetMover()
        cstMover.constraint_file(cst_file)

    # 4. create ScoreFunctions for centroid and fullatom docking
    scorefxn = create_score_function('dna')
    scorefxn.set_weight(core.scoring.fa_elec , 1)    # an "electrostatic" term
    scorefxn.set_weight(core.scoring.atom_pair_constraint, 1.0) #enable our constraints

    #### global docking, a problem solved by the Rosetta DockingProtocol,
    ####    requires interface detection and refinement
    #### as with other protocols, these tasks are split into centroid (interface
    ####    detection) and high-resolution (interface refinement) methods
    #### without a centroid representation, low-resolution DNA-protein
    ####    prediction is not possible and as such, only the high-resolution
    ####    DNA-protein interface refinement is available
    #### WARNING: if you add a perturbation or randomization step, the
    ####    high-resolution stages may fail (see Changing DNA Docking
    ####    Sampling below)
    #### a perturbation step CAN make this a global docking algorithm however
    ####    the rigid-body sampling preceding refinement will require EXTENSIVE
    ####    sampling to produce accurate results and this algorithm spends most
    ####    of its effort in refinement (which may be useless for the predicted
    ####    interface)

    # 5. setup the high resolution (fullatom) docking protocol (DockMCMProtocol)
    # ...as should be obvious by now, Rosetta applications have no central
    #    standardization, the DockingProtocol object can be created and
    #    applied to perform Rosetta docking, many of its options and settings
    #    can be set using the DockingProtocol setter methods
    # as there is currently no centroid representation of DNA in the chemical
    #    database, the low-resolution docking stages are not useful for
    #    DNA docking
    # instead, create an instance of just the high-resolution docking stages
    docking = protocols.docking.DockMCMProtocol()
    docking.set_scorefxn(scorefxn)

    # 6. setup the PyJobDistributor
    jd = PyJobDistributor( job_output , jobs , scorefxn )

    # 7. setup a PyMOL_Observer (optional)
    # the PyMOL_Observer object owns a PyMOL_Mover and monitors pose objects for
    #    structural changes, when changes are detected the new structure is
    #    sent to PyMOL
    # fortunately, this allows investigation of full protocols since
    #    intermediate changes are displayed, it also eliminates the need to
    #    manually apply the PyMOL_Mover during a custom protocol
    # unfortunately, this can make the output difficult to interpret (since you
    #    aren't explicitly telling it when to export) and can significantly slow
    #    down protocols since many structures are output (PyMOL can also slow
    #    down if too many structures are provided and a fast machine may
    #    generate structures too quickly for PyMOL to read, the
    #    "Buffer clean up" message
    if use_pymol:
        pyobs = protocols.moves.AddPyMolObserver(test_pose,True,0.0)

    # 8. perform protein-protein docking
    counter = 0    # for pretty output to PyMOL
    while not jd.job_complete:
        # a. set necessary variables for this trajectory
        # -reset the test pose to original (centroid) structure
        test_pose.assign(pose)
        # -change the pose name, for pretty output to PyMOL
        counter += 1
        test_pose.pdb_info().name(job_output + '_' + str(counter))

        if cst_file != '':
            #apply constraints
            cstMover.apply(test_pose)
        # b. perform docking
        docking.apply(test_pose)

        # c. output the decoy structure:
        # to PyMOL
        test_pose.pdb_info().name(job_output + '_' + str(counter) + '_fa')
        # to a PDB file
        jd.output_decoy(test_pose)
        #pymover.apply(test_pose)

def make_pdbs(pdb_filename="test.pdb",t7="t7_base.pdb",px="PX.pdb",phosphate="6_X"):
    """Using the given files t7 and px and the scissile phosphate, this function
    will generate a new pdb file containing both structures with the scissile
    phosphate of px in the position of the scissile phosphate of the Holliday
    junction in T7"""

    #DNA backbone
    backbone = ['P','OP1','OP2','O5\'','C5\'','C4\'','O4\'','C3\'','C4\'','O4\'','C2\'','C1\'']

    #parse scissile phosphate information
    phosphate_resi, phosphate_chain = phosphate.split("_")
    phosphate_resi = int(phosphate_resi)
    print phosphate_chain, phosphate_resi

    #load in structure files
    parser = Bio.PDB.PDBParser()
    px_structure = parser.get_structure('PX',px)
    t7_structure = parser.get_structure('T7',t7)

    #Perform superimposition
    sup = Bio.PDB.Superimposer()


    fixed = [atom for atom in t7_structure[0]['Y'][8]]
    moving = [atom for atom in px_structure[0][phosphate_chain][phosphate_resi]]

    #Get rid of hydrogens which may or may not exist in either structure and keep only
    #the phosphate backbone
    i = 0
    size = len(fixed)
    while i < size:
        if fixed[i].get_id() not in backbone:
            fixed.pop(i)
            size -= 1
        else:
            i += 1
    i = 0
    size = len(moving)
    while i < size:
        if  moving[i].get_id() not in backbone:
            moving.pop(i)
            size -= 1
        else:
            i += 1

    #If moving starts with the first nucleotide, there will be no phosphate at the beginning
    if phosphate_resi == px_structure[0][phosphate_chain].get_list()[0].get_id()[1]:
        for i in range(3):
            fixed.pop(0)
    #Calculate the superimposition matrices
    try:
        sup.set_atoms(fixed,moving)
    except Exception as e:
        print "Error superimposing PX. Atom number mismatch."
        print 'fixed atoms:',len(fixed),fixed
        print 'moving atoms:',len(moving),moving
        raise e
    #Apply the superimposition to all of PX
    sup.apply(px_structure[0].get_atoms())

    #setup final structure
    output_structure = Bio.PDB.Structure.Structure(pdb_filename)
    output_structure.add(Bio.PDB.Model.Model(0))
    for chain in t7_structure[0]:
        #We don't want to add the DNA from the T7 sequence
        if chain.get_id() not in ['Y','Z']:
            chain.detach_parent()
            output_structure[0].add(chain)
    for chain in px_structure[0]:
        chain.detach_parent()
        output_structure[0].add(chain)

    #and output final structure
    io = Bio.PDB.PDBIO()
    io.set_structure(output_structure)
    io.save(pdb_filename)

def rescore(pdb_filename,cst_file):
    '''helper function for rescoring PDBs.  Only meant to be used when imported
    and run from commandline interpreter'''
    pose = Pose()
    if PR4:
        pose_from_file(pose, pdb_filename)
    else:
        pose_from_pdb(pose,pdb_filename)

    cstMover = protocols.simple_moves.ConstraintSetMover()
    cstMover.constraint_file(cst_file)

    # 4. create ScoreFunctions for centroid and fullatom docking
    scorefxn = create_score_function('dna')
    scorefxn.set_weight(core.scoring.fa_elec , 1)    # an "electrostatic" term
    scorefxn.set_weight(core.scoring.atom_pair_constraint, 1.0) #enable our constraints
    cstMover.apply(pose)
    scorefxn.show(pose)
    pose.dump_scored_pdb(pdb_filename,scorefxn)

def create_constraints(cst_filename, phosphate):
    """Generates the constraints file to be used for the docking simulation based
    on the selected scissile phosphate"""
    phosphate = ''.join(phosphate.split("_"))
    #entries in constraint list are tuples in the following format:
    #(atomname, residuenumber, px_atom, dist, dist_cutoff)
    constraint_list = [('CG','55A','P',5.37,0.5),('CD','20B','P',5.53,0.5),('CA','150C','OP1',2.43,0.25),('CA','30C','OP1',2.15,0.25),
        ('CA','30C','P',2.95,0.25),('CA','150C','P',3.38,0.25)]
    cst_out = open(cst_filename,'w+')
    for cst in constraint_list:
        cst_out.write('AtomPair %s %s %s %s HARMONIC %f %f\n'%(cst[0],cst[1],cst[2],phosphate,cst[3],cst[4]))
    cst_out.close()

def sort_fasc_file(fasc_file):
    '''Takes the .fasc file output by sample_dna_interface and rearranges it to be
    in ascending order based on total energy score'''
    first = ''
    with open(fasc_file) as fin:
        lines = [line.strip().split() for line in fin]
        fin.close()
        first = lines[0]
        lines.pop(0)
    output = sorted(lines,key=lambda tup: float(tup[3]))
    fout = open(fasc_file,'w+')
    fout.write(" ".join(first))
    fout.write('\n')
    for line in output:
        fout.write(" ".join(line))
        fout.write('\n')
    fout.close()

def minimize(pdb_in,pdb_out,cst_file):
    '''minimizes the pdb file pdb_in using the same scorefunction and constraints file used for docking
    outputs minimized pdb as pdb_out with score information'''
    pose = Pose()
    if PR4:
        pose_from_file(pose,pdb_in)
    else:
        pose_from_pdb(pose,pdb_in)
    #scorefunction setup
    sfxn = create_score_function("dna")
    sfxn.set_weight(core.scoring.fa_elec,1.0)
    sfxn.set_weight(core.scoring.atom_pair_constraint,1.0)
    #cstMover setup
    cstMover = protocols.simple_moves.ConstraintSetMover()
    cstMover.constraint_file(cst_file)
    cstMover.apply(pose)
    #movemap setup
    movemap = core.kinematics.MoveMap()
    movemap.set_bb_true_range(1,133)
    movemap.set_chi_true_range(1,133)
    for i in range(1,8):
        if i < 5:
            movemap.set_jump(i,True)
        else:
            movemap.set_jump(i,False)
    #minMover setup
    minMover = protocols.simple_moves.MinMover()
    minMover.score_function(sfxn)
    minMover.movemap(movemap)
    #print score before min
    print "score before minimization"
    sfxn.show(pose)
    #do the minimization
    minMover.apply(pose)
    #print score
    sfxn.show(pose)
    #dump pdb
    pose.dump_scored_pdb(pdb_out,sfxn)


def main():
    parser = optparse.OptionParser()
    #source of T7 file
    parser.add_option('--t7',dest='t7',default='t7_base.pdb',help='PDB file specifying the coordinates of the single-headed T7 model')
    #source of PX file
    parser.add_option('--px',dest='px',default='PX.pdb',help='PDB file specifying the coordinates of the PX DNA model')
    #Scissile phosphate
    parser.add_option('--phosphate',dest='phosphate',default='6_X',help='String specifying the scissile phosphate of PX to model.  This should be in the format nucleotide#_chainID (e.g. 6_X for 6th nucleotide of chain X)')
    #pymol flag
    parser.add_option('--pymol',dest='use_pymol',action='store_true',default=False,help='use this flag if you want to view live output in pymol')
    #jobs
    parser.add_option('--jobs',type='int',dest='jobs',default=1,help='the number of jobs (trajectories) to perform')
    #job_output
    parser.add_option('--job_output',dest='job_output',default='dna_output',help='the name preceding all output, output PDB files and .fasc')

    (options,args) = parser.parse_args()
    pdb_filename = "%s.pdb"%(options.job_output)
    jobs = options.jobs
    use_pymol = options.use_pymol
    phosphate = options.phosphate.upper()
    # try:
    #     os.chdir(options.job_output)
    # except:
    #     os.mkdir(options.job_output)
    #     os.chdir(options.job_output)
    #     shutil.copy("%s/../%s"%(os.get_cwd(),options.t7),"%s/"%(os.get_cwd()))
    #     shutil.copy("%s/../%s"%(os.get_cwd(), options.px),"%s/"%(os.get_cwd()))

    cst_filename = "%s.cst"%(options.job_output)
    make_pdbs(pdb_filename,options.t7,options.px,phosphate)
    create_constraints(cst_filename,phosphate)
    sample_dna_interface(pdb_filename, 'ABC_WXYZ', jobs, options.job_output, cst_filename, use_pymol )
    sort_fasc_file("%s.fasc"%(options.job_output))

if __name__ == '__main__':
    main()

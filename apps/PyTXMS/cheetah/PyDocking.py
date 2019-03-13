#!usr/bin/env python

from __future__ import print_function
from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *

def PyDocking(pose, partners, translation, rotation, jobs, job_output, low_only, high_only):

    print ('Docking partners selected as', partners)
    dock_jump = 1
    # movable_jmp = Vector1([partners.find('_')])
    # print ("movable_jmp is: ", movable_jmp)

    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')

    recover_sidechains = protocols.simple_moves.ReturnSidechainMover(pose)

    to_centroid.apply(pose)

    test_pose = Pose()
    test_pose.assign(pose)

    # scorefxn_low = create_score_function('interchain_cen')
    # scorefxn_high = create_score_function('docking')
    # scorefxn_high_min = create_score_function('docking', 'docking_min')

    randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump,
        partner_upstream)
    randomize_downstream = RigidBodyRandomizeMover(pose, dock_jump,
        partner_downstream)

    dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
    spin = RigidBodySpinMover(dock_jump)

    slide_into_contact = protocols.docking.DockingSlideIntoContact(dock_jump)


    movemap = MoveMap()
    movemap.set_jump(dock_jump, True)

    minmover = protocols.minimization_packing.MinMover()
    minmover.movemap(movemap)
    # minmover.score_function(scorefxn_high_min)

    if low_only:
        perturb = protocols.moves.SequenceMover()
        perturb.add_mover(randomize_upstream)
        perturb.add_mover(randomize_downstream)
        perturb.add_mover(dock_pert)
        perturb.add_mover(spin)
        perturb.add_mover(slide_into_contact)
        perturb.add_mover(to_fullatom)
        perturb.add_mover(recover_sidechains)
        perturb.add_mover(minmover)
        print('low resolution docking ...\n')
        dock_prot = protocols.docking.DockingProtocol()
        dock_prot.set_low_res_protocol_only(True)
        dock_prot.set_movable_jumps(Vector1([1]))
        # dock_prot.set_lowres_scorefxn(scorefxn_low)
        # jd = PyJobDistributor(job_output, jobs, scorefxn_low)
        # temp_pose = Pose()
        # temp_pose.assign(pose)
        # jd.native_pose = temp_pose

    elif high_only:
        perturb = protocols.moves.SequenceMover()
        # perturb.add_mover(randomize_upstream)
        # perturb.add_mover(randomize_downstream)
        perturb.add_mover(dock_pert)
        # perturb.add_mover(spin)
        perturb.add_mover(slide_into_contact)
        perturb.add_mover(to_fullatom)
        perturb.add_mover(recover_sidechains)
        perturb.add_mover(minmover)
        print('high resolution docking ...\n')
        dock_prot = protocols.docking.DockingProtocol()
        dock_prot.set_docking_local_refine(True)
        dock_prot.set_movable_jumps(Vector1([1]))
        # dock_prot.set_highres_scorefxn(scorefxn_high_min)
        # jd = PyJobDistributor(job_output, jobs, scorefxn_high)
        # temp_pose = Pose()
        # temp_pose.assign(pose)
        # to_fullatom.apply(temp_pose)
        # recover_sidechains.apply(temp_pose)
        # jd.native_pose = temp_pose

    else:
        perturb = protocols.moves.SequenceMover()
        perturb.add_mover(randomize_upstream)
        perturb.add_mover(randomize_downstream)
        perturb.add_mover(dock_pert)
        perturb.add_mover(spin)
        perturb.add_mover(slide_into_contact)
        perturb.add_mover(to_fullatom)
        perturb.add_mover(recover_sidechains)
        perturb.add_mover(minmover)
        print('default docking protocol...\n')
        dock_prot = protocols.docking.DockingProtocol()
        dock_prot.set_movable_jumps(Vector1([1]))
        # dock_prot.set_lowres_scorefxn(scorefxn_low)
        # dock_prot.set_highres_scorefxn(scorefxn_high_min)
        # jd = PyJobDistributor(job_output, jobs, scorefxn_high)
        # temp_pose = Pose()
        # temp_pose.assign(pose)
        # to_fullatom.apply(temp_pose)
        # recover_sidechains.apply(temp_pose)
        # jd.native_pose = temp_pose



    # perform protein-protein docking
    # counter = 0
    # while not jd.job_complete:

    try:
        test_pose.assign(pose)
        # counter += 1
        # test_pose.pdb_info().name(job_output + 'model_' + str(counter))
        perturb.apply(test_pose)

        if (low_only):
            dock_prot.apply(test_pose)
        elif (high_only):
            dock_prot.apply(test_pose)
            to_fullatom.apply(test_pose)
        else:
            dock_prot.apply(test_pose)
            to_fullatom.apply(test_pose)

        # test_pose.pdb_info().name(job_output + 'model_' + str( counter ) + '_fa')

        return test_pose
        # jd.output_decoy(test_pose)
    except RuntimeError:
        print('ERROR:Error occurred in calculations!')



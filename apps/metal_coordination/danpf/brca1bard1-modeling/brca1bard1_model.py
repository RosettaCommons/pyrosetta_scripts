#!/usr/bin/env python
"""
usage:
    python brca1bard1_model.py 1

cite:

BRCA1/BARD1 site-specific ubiquitylation of nucleosomal H2A is directed by BARD1.
Sam Witus, Anika Burrell, Daniel Farrell, Jianming Khang, Meiling Wang, Jesse Hansen, Alex Pravat, Lisa Tuttle, Peter Brzovic, Champak Chatterjee, Weixing Zhao, Frank DiMaio, Justin Kollman, and Rachel Klevit

"""

import sys
import os
import pyrosetta
import pathlib
import itertools
from typing import List

from pyrosetta import (
    rosetta,
    MoveMap,
    init,
    create_score_function,
    pose_from_pdb,
    standard_task_factory,
)

from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
import copy


def remove_h(pose, resn, atomn):
    res = pose.residue(resn)
    current_residue_type_basename = rosetta.core.chemical.residue_type_base_name(res.type())
    current_residue_type_patches_name = rosetta.core.chemical.residue_type_all_patches_name(res.type())
    res_atom = res.atom_name(atomn)
    res_patchname = f"MP-{res_atom.strip()}-pruneH"
    if res_patchname in current_residue_type_patches_name:
        patches = current_residue_type_patches_name.split(":")
        res_type_mod_name = copy.copy(current_residue_type_basename)
        found = False
        for patch in patches:
            if not patch:
                continue
            if patch == res_patchname and found is False:
                res_type_mod_name = f"{res_type_mod_name}:{patch}"
                found = True
            elif patch != res_patchname:
                res_type_mod_name = f"{res_type_mod_name}:{patch}"

    else:
        res_type_mod_name = f"{current_residue_type_basename}:{res_patchname}{current_residue_type_patches_name}"
    restype_for_pose = rosetta.core.pose.get_restype_for_pose(pose, res_type_mod_name, res.type().mode())
    new_res = rosetta.core.conformation.Residue(restype_for_pose, True)
    old_res = res.clone()
    pose.replace_residue(resn, new_res, True)
    for x in range(1, old_res.natoms() + 1):
        if new_res.has(old_res.atom_name(x)):
            pose.set_xyz(
                rosetta.core.id.AtomID(new_res.atom_index(old_res.atom_name(x)), resn),
                old_res.xyz(old_res.atom_name(x)),
            )


def get_close_residues(pose, center_xyz, residues_to_skip: List[int]):
    close_residues = []
    for resnum in range(1, pose.size() + 1):
        if resnum in residues_to_skip:
            continue
        res = pose.residue(resnum)
        for atomnum in range(1, res.nheavyatoms() + 1):
            if atomnum == 5:  # CB
                distance = res.atom(atomnum).xyz().distance(center_xyz)
                if distance < 7:
                    close_residues.append((resnum, distance))

    close_residues.sort(key=lambda x: x[1])
    return close_residues


def setup_zn_coordination_center(pose, zn_residue: int, close_residues_w_distance):
    zn_xyz = pose.residue(zn_residue).atom(1).xyz()
    zn_atom_id = rosetta.core.id.AtomID(1, zn_residue)
    count = 0
    constraints = []

    tet_atom_ids = []
    # Make sure you keep the correct atom_id after deleting some hydrogens!
    for current_residue_num, cb_distance in close_residues_w_distance:
        if pose.residue(current_residue_num).name3() == "CYS":
            atom_index = pose.residue(current_residue_num).atom_index(" SG ")
            remove_h(pose, current_residue_num, atom_index)
            fn = rosetta.core.scoring.func.HarmonicFunc(2.34, 0.1)
            this_atom_id = rosetta.core.id.AtomID(
                pose.residue(current_residue_num).atom_index(" SG "), current_residue_num
            )
            cst = rosetta.core.scoring.constraints.AtomPairConstraint(
                rosetta.core.id.AtomID(1, zn_residue), this_atom_id, fn
            )
            constraints.append(cst)
            cys_angle_fn = rosetta.core.scoring.func.CircularHarmonicFunc(1.90241, 0.1)
            angle_cst = rosetta.core.scoring.constraints.AngleConstraint(
                zn_atom_id, this_atom_id, rosetta.core.id.AtomID(5, current_residue_num), cys_angle_fn
            )
            constraints.append(angle_cst)
            tet_atom_ids.append(this_atom_id)
            count += 1
        if pose.residue(current_residue_num).name3() == "HIS":
            # Just take closer taut
            nd1_atom_index = pose.residue(current_residue_num).atom_index(" ND1")
            nd1_atom_distance = pose.residue(current_residue_num).atom(nd1_atom_index).xyz().distance(zn_xyz)
            ne2_atom_index = pose.residue(current_residue_num).atom_index(" NE2")
            ne2_atom_distance = pose.residue(current_residue_num).atom(ne2_atom_index).xyz().distance(zn_xyz)
            if nd1_atom_distance < ne2_atom_distance:
                remove_h(pose, current_residue_num, nd1_atom_index)
                atom_index = nd1_atom_index
                atom_index = pose.residue(current_residue_num).atom_index(" ND1")
                print(current_residue_num, cb_distance, "found his", "ND1")
            else:
                remove_h(pose, current_residue_num, ne2_atom_index)
                atom_index = pose.residue(current_residue_num).atom_index(" NE2")
                print(current_residue_num, cb_distance, "found his", "NE2")
            current_residue = pose.residue(current_residue_num)  # get after removing h

            this_atom_id = rosetta.core.id.AtomID(atom_index, current_residue_num)
            fn = rosetta.core.scoring.func.HarmonicFunc(2.057, 0.1)
            cst = rosetta.core.scoring.constraints.AtomPairConstraint(
                rosetta.core.id.AtomID(1, zn_residue), this_atom_id, fn
            )
            constraints.append(cst)
            zero_dihedral_fn = rosetta.core.scoring.func.CircularHarmonicFunc(0, 0.1)
            angle_fn = rosetta.core.scoring.func.CircularHarmonicFunc(2.0944, 0.1)
            if atom_index == nd1_atom_index:
                dihed_cst = rosetta.core.scoring.constraints.DihedralConstraint(
                    rosetta.core.id.AtomID(current_residue.atom_index(" CB "), current_residue_num),
                    rosetta.core.id.AtomID(current_residue.atom_index(" CG "), current_residue_num),
                    this_atom_id,
                    zn_atom_id,
                    zero_dihedral_fn,
                )
                angle_cst = rosetta.core.scoring.constraints.AngleConstraint(
                    rosetta.core.id.AtomID(current_residue.atom_index(" CE1"), current_residue_num),
                    this_atom_id,
                    zn_atom_id,
                    angle_fn
                )
                # Angle CE1 1060 ND1 1060 ZN 1263 CIRCULARHARMONIC 2.0944 0.1
            else:
                dihed_cst = rosetta.core.scoring.constraints.DihedralConstraint(
                    rosetta.core.id.AtomID(current_residue.atom_index(" HD2"), current_residue_num),
                    rosetta.core.id.AtomID(current_residue.atom_index(" CD2"), current_residue_num),
                    this_atom_id,
                    zn_atom_id,
                    zero_dihedral_fn,
                )
                angle_cst = rosetta.core.scoring.constraints.AngleConstraint(
                    rosetta.core.id.AtomID(current_residue.atom_index(" CD2"), current_residue_num),
                    this_atom_id,
                    zn_atom_id,
                    angle_fn
                )
                # Angle CE1 801 NE2 801 ZN 1261 CIRCULARHARMONIC 2.0944 0.1
            constraints.append(angle_cst)
            constraints.append(dihed_cst)
            tet_atom_ids.append(this_atom_id)
            count += 1
        if count == 4:  # Tetrahedral
            break
    # setup tetrahedral constraints
    tet_angle_fn = rosetta.core.scoring.func.CircularHarmonicFunc(1.9111355, 0.1)
    tet_size = 0
    for combination in itertools.combinations(tet_atom_ids, 2):
        cst = rosetta.core.scoring.constraints.AngleConstraint(
            combination[0], rosetta.core.id.AtomID(1, zn_residue), combination[1], tet_angle_fn
        )
        constraints.append(cst)
        tet_size += 1
    for constraint in constraints:
        if constraint.type() == "AtomPair":
            ostream = pyrosetta.rosetta.std.ostringstream()
            cst.show_def(ostream, pose)
            print(ostream.str())
            continue
        print(constraint, end="")
    return constraints


def setup_zn_coordination_centers(pose):
    zn_residues = []
    for resnum in range(1, pose.size() + 1):
        if pose.residue(resnum).name3() == " ZN":
            zn_residues.append(resnum)

    cst_set = rosetta.core.scoring.constraints.ConstraintSet()
    for zn_residue in zn_residues:
        constraints = setup_zn_coordination_center(
            pose, zn_residue, get_close_residues(pose, pose.residue(zn_residue).atom(1).xyz(), [zn_residue])
        )
        for constraint in constraints:
            cst_set.add_constraint(constraint)
            print(constraint)

    less_than_3_fn = rosetta.core.scoring.func.FlatHarmonicFunc(0, 0.1, 2.9)
    # Salt bridge 01
    nz_71 = rosetta.core.id.AtomID(
        pose.residue(71).atom_index(" NZ "), 71
    )
    oe1_409 = rosetta.core.id.AtomID(
        pose.residue(409).atom_index(" OE1"), 409
    )
    cst = rosetta.core.scoring.constraints.AtomPairConstraint(
        nz_71, oe1_409, less_than_3_fn
    )
    cst_set.add_constraint(cst)

    # Salt bridge 02
    less_than_3_fn2 = rosetta.core.scoring.func.FlatHarmonicFunc(0, 0.1, 3.0)
    nh1_722 = rosetta.core.id.AtomID(
        pose.residue(437).atom_index(" OE2"), 437
    )
    oe2_4372 = rosetta.core.id.AtomID(
        pose.residue(538).atom_index(" N  "), 538
    )
    cst2 = rosetta.core.scoring.constraints.AtomPairConstraint(
        nh1_722, oe2_4372, less_than_3_fn2
    )
    cst_set.add_constraint(cst2)
    # Salt bridge 03
    less_than_3_fn2 = rosetta.core.scoring.func.FlatHarmonicFunc(0, 0.1, 2.9)
    oe2_437 = rosetta.core.id.AtomID(
        pose.residue(437).atom_index(" OE2"), 437
    )
    nh1_72 = rosetta.core.id.AtomID(
        pose.residue(72).atom_index(" NH1"), 72
    )
    cst3 = rosetta.core.scoring.constraints.AtomPairConstraint(
        oe2_437, nh1_72, less_than_3_fn2
    )
    cst_set.add_constraint(cst3)

    # Salt bridge 04
    less_than_32_fn2 = rosetta.core.scoring.func.FlatHarmonicFunc(0, 0.1, 3.0)
    nz_39 = rosetta.core.id.AtomID(
        pose.residue(39).atom_index(" NZ "), 39
    )
    oe1_538 = rosetta.core.id.AtomID(
        pose.residue(538).atom_index(" OE2"), 538
    )
    cst4 = rosetta.core.scoring.constraints.AtomPairConstraint(
        nz_39, oe1_538, less_than_32_fn2
    )
    cst_set.add_constraint(cst4)

    # Setting angles properly
    angle_170 = rosetta.core.scoring.func.CircularHarmonicFunc(2.96706, 0.1)
    oe2_437 = rosetta.core.id.AtomID(
        pose.residue(437).atom_index(" OE2"), 437
    )
    N_538 = rosetta.core.id.AtomID(
        pose.residue(538).atom_index(" N  "), 538
    )
    H_538 = rosetta.core.id.AtomID(
        pose.residue(538).atom_index(" H  "), 538
    )
    cst_ang = rosetta.core.scoring.constraints.AngleConstraint(
        N_538, H_538, oe2_437, angle_170
    )
    cst_set.add_constraint(cst_ang)

    cstm = rosetta.protocols.constraint_movers.ConstraintSetMover()
    cstm.constraint_set(cst_set)
    cstm.apply(pose)


def atom_ids_from_names(names, pose, n_by_olc):
    split = names.split("-")
    atom_ids = []
    print(names)
    for unsplit_id in split:
        print(unsplit_id)
        olc, atom_name = unsplit_id.split("_")
        atom_ids.append(
            rosetta.core.id.AtomID(
                pose.residue(n_by_olc[olc]).atom_index(f" {atom_name} "), n_by_olc[olc]))
    return atom_ids


def main():
    out_folder = "output_001"
    pathlib.Path(out_folder).mkdir(exist_ok=True)
    pose_dump_fn = f"{out_folder}/minrun_{int(sys.argv[1]):03}"
    if os.path.isfile(pose_dump_fn):
        sys.exit(-1)
    todel = []
    for key in os.environ.keys():
        if "ROSETTA" in key:
            todel.append(key)

    for key in todel:
        del os.environ[key]

    init_flags = (
        "-default_max_cycles 200"
        " -missing_density_to_jump"
        " -edensity:mapreso 2.5"
        " -edensity:mapfile ./denmod_resamp.mrc"
        " -detect_disulf false"
        " -out:pdb"
        " -chemical:set_atomic_charge fa_standard:ZN:ZN:0.0"
        " -relax:ramady true"
        " -relax:ramady_force true"
        " -newdna true"
        " -exclude_dna_dna false"
        " -relax:dna_move true"
        " -beta"
        " -beta_cart"
        " -ex1 "
    )
    print(init_flags)

    init(init_flags)

    cart_sfxn = create_score_function("beta_cart")
    cart_sfxn.set_weight(rosetta.core.scoring.elec_dens_fast, 30)
    cart_sfxn.set_weight(rosetta.core.scoring.cart_bonded, 1.5)
    # Set this high because of close N-Zn in his
    cart_sfxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 5)
    cart_sfxn.set_weight(rosetta.core.scoring.angle_constraint, 5)
    cart_sfxn.set_weight(rosetta.core.scoring.dihedral_constraint, 10)
    cart_sfxn.set_weight(rosetta.core.scoring.coordinate_constraint, 5)

    c_mmap = MoveMap()
    c_mmap.set_bb(False)
    c_mmap.set_chi(False)
    c_mmap.set_jump(False)
    c_mmap.set(rosetta.core.id.D, True)

    min_mover_d = MinMover(c_mmap, cart_sfxn, "lbfgs_armijo_nonmonotone", 0.0001, True)
    min_mover_d.max_iter(200)
    min_mover_d.cartesian(True)

    pose = pose_from_pdb("input_model.pdb")
    rosetta.core.pose.addVirtualResAsRoot(pose)
    setup_zn_coordination_centers(pose)

    lr_cart = rosetta.protocols.relax.LocalRelax()
    lr_cart.set_sfxn(cart_sfxn)
    lr_cart.set_max_iter(200)
    lr_cart.set_ncyc(4)
    lr_cart.cartesian(True)

    # initial packing
    tf = standard_task_factory()
    task = tf.create_task_and_apply_taskoperations(pose)
    task.restrict_to_repacking()
    task.temporarily_set_pack_residue(591, True)
    task.temporarily_set_pack_residue(437, True)
    task.temporarily_set_pack_residue(595, True)
    task.temporarily_set_pack_residue(706, True)
    task.temporarily_set_pack_residue(39, True)
    task.temporarily_set_pack_residue(73, True)
    task.temporarily_set_pack_residue(538, True)
    pack_mover = PackRotamersMover(cart_sfxn, task)
    pack_mover.apply(pose)

    min_mover_d.apply(pose)
    pose.dump_pdb(f"{pose_dump_fn}_packmm_t01.pdb")
    for x in range(10):
        lr_cart.apply(pose)
        pose.dump_pdb(f"{pose_dump_fn}_lr_t{x:02}.pdb")
        pack_mover.apply(pose)
        pose.dump_pdb(f"{pose_dump_fn}_pm_t{x:02}.pdb")


if __name__ == "__main__":
    main()

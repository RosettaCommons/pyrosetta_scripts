#!usr/bin/env python

import pyrosetta

from collections import defaultdict


pyrosetta.init()

pose = pyrosetta.pose_from_file("test_pose.pdb")

scorefxn = pyrosetta.get_fa_scorefxn()
scorefxn(pose)

sela = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
selb = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")

neigha = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(sela, 10, False)
neighb = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(selb, 10, False)

suba = neigha.apply(pose)
subb = neighb.apply(pose)


# this will return all currently active TwoBodyEnergies
# this thing returns two AnalyticEtableEnergy which is correct.
# they get initialized with different flags for intra and inter
def find_all_2body_methods(scorefxn):
    nonzero = scorefxn.get_nonzero_weighted_scoretypes()

    to_find = {}
    for scoretype in nonzero:
        to_find[scoretype] = False

    sm = pyrosetta.rosetta.core.scoring.ScoringManager.get_instance()

    found = []

    for scoretype in to_find:
        if (to_find[scoretype]):
            continue

        method = sm.energy_method( scoretype, scorefxn.energy_method_options())
        if (isinstance(method, pyrosetta.rosetta.core.scoring.methods.TwoBodyEnergy)):
            score_types = method.score_types()
            if (scoretype in score_types):
                for this_score_type in method.score_types():
                    if (this_score_type in to_find):
                        to_find[this_score_type] = True  
                found.append(method)

    return found





#this thing returns a dict that takes seqpos ints to access
# elements
# dict[seqposa][seqposb] = emap
# remember to multiply output by scoreterm weights
def find_res_by_res_scores(pose, scorefxn, subset_a, subset_b):

    res_a = pyrosetta.rosetta.core.select.get_residues_from_subset(subset_a)
    res_b = pyrosetta.rosetta.core.select.get_residues_from_subset(subset_b)

    # Make sure hydrogen bonds show up correctly
    # We don't actually have to clone and call set_energy_method_options
    #  But technically energy_method_options() returns a const & and
    #  this is "more correct"
    options = scorefxn.energy_method_options().clone()
    hb_opts = options.hbond_options()
    hb_opts.decompose_bb_hb_into_pair_energies(True) 
    options.hbond_options(hb_opts)
    scorefxn.set_energy_method_options(options)


    twobody_methods = find_all_2body_methods(scorefxn)

    output_dict = defaultdict( lambda : {}, {})

    for i in res_a:
        resi = pose.residue(i)
        for j in res_b:
            resj = pose.residue(j)

            emap = pyrosetta.rosetta.core.scoring.EMapVector()

            for method in twobody_methods:
                method.residue_pair_energy(resi, resj, pose, scorefxn, emap)

            output_dict[i][j] = emap

    return dict(output_dict)    # kill default dict to prevent accidents later



find_res_by_res_scores(pose, scorefxn, suba, subb)
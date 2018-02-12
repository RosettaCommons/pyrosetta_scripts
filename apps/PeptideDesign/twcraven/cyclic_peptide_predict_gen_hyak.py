from pyrosetta import *
import numpy as np
import random, math
from rosetta.protocols.cyclic_peptide import DeclareBond
from rosetta.core.scoring.func import HarmonicFunc

init('-database /suppscr/baker/twcraven/database -override_rsd_type_limit -ex1 -ex2 -ignore_unrecognized_res -output_virtual -mute all -bbdep_omega false -corrections::beta_nov16')

create_residue = rosetta.core.conformation.ResidueFactory.create_residue
chm = rosetta.core.chemical.ChemicalManager.get_instance()
sm = rosetta.core.scoring.ScoringManager.get_instance()

rama = sm.get_RamaPrePro()
#help(rama)


def dist(x, y):
	return math.sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2] - y[2])**2)

rts = chm.residue_type_set( 'fa_standard' )
sfxn = create_score_function('beta_cst')
sfxn2 = create_score_function('beta_nov16')
sfxn3 = create_score_function('mm_std')
#help(sfxn)

movemap = MoveMap()
movemap.set_chi(True)
movemap.set_bb(True)
minmover = rosetta.protocols.simple_moves.MinMover(movemap, sfxn, "linmin", 0.001, True )


for kk in range(0, int(sys.argv[3]))
	res = rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY_p:N_Methylation_noHs'))
	pose = Pose()
	res_length = int(sys.argv[1])
	for i in range(0, res_length):
		if i == 1:
			if sys.argv[2] == 'SAR':
				pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY_p:N_Methylation_noHs')), True)
			elif sys.argv[2] == 'beta':
				pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('B3A')), True)
			elif sys.argv[2] == 'alpha':
				pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA')), True)
			elif sys.argv[2] == 'oligourea':
				pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('BB0')), True)
		elif i == 2 and sys.argv[2] == 'oligourea':
			pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY')), True)


		else:
			pose.append_residue_by_bond(res, True)
	

	if sys.argv[2] == 'oligourea':

		for i in range(1, res_length):
			torsions = pyrosetta.rosetta.utility.vector1_double()
			
			if i != 2:
				rama.random_mainchain_torsions(pose.conformation(), pose.residue(i).type(), pose.residue(i+1).type(), torsions)
			else:
				torsions = (random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360))
			
			if i != 3 and i != 1:
				ri = np.random.randint(2, size=1)
				if ri == 1:
					pose.set_omega(i, 180.0)
				else:
					pose.set_omega(i, 0.0)
			elif i == 1:
				pose.conformation().set_torsion_angle(AtomID(pose.residue(2).atom_index('CA'), 2), AtomID(pose.residue(2).atom_index('N'), 2), AtomID(pose.residue(1).atom_index('C'), 1), AtomID(pose.residue(1).atom_index('CA'), 1),  180*0.0174532925)
			else:
				pose.conformation().set_torsion_angle(AtomID(pose.residue(3).atom_index('CA'), 3), AtomID(pose.residue(3).atom_index('N'), 3), AtomID(pose.residue(2).atom_index('C'), 2), AtomID(pose.residue(2).atom_index('NU'), 2),  180*0.0174532925)

			if i != 2:
				for j in range(1, len(torsions)):
					pose.set_torsion( rosetta.core.id.TorsionID(i, rosetta.core.id.BB, j), torsions[j] )
			else:
				pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('C'), 1), AtomID(pose.residue(2).atom_index('N'), 2), AtomID(pose.residue(2).atom_index('CA'), 2), AtomID(pose.residue(2).atom_index('CM'), 2), torsions[0]*0.0174532925)
				pose.conformation().set_torsion_angle(AtomID(pose.residue(2).atom_index('N'), 2), AtomID(pose.residue(2).atom_index('CA'), 2), AtomID(pose.residue(2).atom_index('CM'), 2), AtomID(pose.residue(2).atom_index('NU'), 2),  torsions[1]*0.0174532925)
				pose.conformation().set_torsion_angle(AtomID(pose.residue(2).atom_index('CA'), 2), AtomID(pose.residue(2).atom_index('CM'), 2), AtomID(pose.residue(2).atom_index('NU'), 2), AtomID(pose.residue(2).atom_index('C'), 2),  torsions[2]*0.0174532925)

	
	else:
		for i in range(1, res_length):
			torsions = pyrosetta.rosetta.utility.vector1_double()

			if i != 2:
				rama.random_mainchain_torsions(pose.conformation(), pose.residue(i).type(), pose.residue(i+1).type(), torsions)
			elif i == 2 and sys.argv[2] == 'beta':
				torsions = (random.uniform(0, 360), random.uniform(0, 360), random.uniform(0, 360))
			elif (i == 2 and sys.argv[2] == 'alpha') or (i == 2 and sys.argv[2] == 'SAR'):
				rama.random_mainchain_torsions(pose.conformation(), pose.residue(i).type(), pose.residue(i+1).type(), torsions)
			
			ri = np.random.randint(2, size=1)
			if ri == 1:
				pose.set_omega(i, 180.0)
			else:
				pose.set_omega(i, 0.0)
				
			for j in range(1, len(torsions)):
				pose.set_torsion( rosetta.core.id.TorsionID(i, rosetta.core.id.BB, j), torsions[j] )
			
	pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('CN'), 1), AtomID(pose.residue(1).atom_index('N'), 1), AtomID(pose.residue(1).atom_index('CA'), 1), AtomID(pose.residue(1).atom_index('C'), 1), random.uniform(0, 360)*0.0174532925)
	pose.conformation().set_torsion_angle(AtomID(pose.residue(res_length).atom_index('O'), res_length), AtomID(pose.residue(res_length).atom_index('C'), res_length), AtomID(pose.residue(res_length).atom_index('CA'), res_length), AtomID(pose.residue(res_length).atom_index('N'), res_length), random.uniform(-125.0, 125.0)*0.0174532925)

	pose.conformation().declare_chemical_bond(1, 'N', res_length, 'C')

	genkic = pyrosetta.rosetta.protocols.generalized_kinematic_closure.GeneralizedKIC()

	cyclization_point_start = res_length 
	cyclization_point_end = 1
	anchor_res_min = 2
	anchor_res_max = res_length -1

	anchor_res = np.random.randint(anchor_res_min, high= anchor_res_max, size = 1)[0]
	first_loop_res =  anchor_res + 1
	last_loop_res = anchor_res - 1 

	middle_loop_res = np.random.randint(cyclization_point_end, high=cyclization_point_start-3, size = 1)[0]
	print anchor_res, middle_loop_res
	if middle_loop_res == last_loop_res:
		middle_loop_res += 3
	elif middle_loop_res == anchor_res:
		middle_loop_res +=2
	elif middle_loop_res == first_loop_res:
		middle_loop_res +=1
	if middle_loop_res > res_length :
		middle_loop_res = middle_loop_res - res_length 

	for i in range(first_loop_res, cyclization_point_start+1):
		genkic.add_loop_residue(i)
	for i in range (cyclization_point_end, last_loop_res+1):
		genkic.add_loop_residue(i)

	genkic.add_perturber('randomize_backbone_by_rama_prepro')
	genkic.add_perturber('sample_cis_peptide_bond')
	genkic.add_value_to_perturber_value_list( 2, 0.5)

	i = anchor_res + 1
	while (i!=anchor_res):
		if i > cyclization_point_start:
			i = cyclization_point_end
		elif i == anchor_res:
			i = i#do nothing
		else:
			genkic.add_residue_to_perturber_residue_list(1, i)
			if i != res_length:
				if pose.residue(i+1).has_property('N_METHYLATED'):
					genkic.add_residue_to_perturber_residue_list(2, i)
		i += 1

	genkic.add_filter('loop_bump_check')

	ri = np.random.randint(2, size=1)
	if ri == 1:
		genkic.close_bond(res_length, 'C', 1, 'N', 0, '', 0, '', 1.328, 116.206, 121.702, 180, False, False)

	else:
		genkic.close_bond(res_length, 'C', 1, 'N', 0, '', 0, '', 1.328, 116.206, 121.702, 0, False, False)

	genkic.set_pivot_atoms( first_loop_res, 'CA', middle_loop_res, 'CA', last_loop_res, 'CA' );
	genkic.set_closure_attempts(2000)
	genkic.set_selector_scorefunction(sfxn3)
	genkic.set_selector_type('lowest_energy_selector')
	genkic.apply(pose)

	for i in range(0, res_length):
		if pose.residue(i+1).has_property('N_METHYLATED'):
			pose.replace_residue(i+1, rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY_p:N_Methylation')), True)

	dhf = HarmonicFunc( 180.0*0.0174533 , 0.1 )
	dhc = rosetta.core.scoring.constraints.DihedralConstraint( AtomID(pose.residue(1).atom_index('N'), 1) , AtomID(pose.residue(res_length).atom_index('CA'), res_length) , AtomID(pose.residue(res_length).atom_index('C'), res_length) , AtomID(pose.residue(res_length).atom_index('O'), res_length) , dhf )
	pose.add_constraint( dhc )

	dhf = HarmonicFunc( 180.0*0.0174533 , 0.1 )
	dhc = rosetta.core.scoring.constraints.DihedralConstraint( AtomID(pose.residue(1).atom_index('CA'), 1) , AtomID(pose.residue(res_length).atom_index('C'), res_length) , AtomID(pose.residue(1).atom_index('N'), 1) , AtomID(pose.residue(1).atom_index('CN'), 1) , dhf )
	pose.add_constraint( dhc )


	min_e = 10**10; best = 0.0
	for i in range(0, 360):
		pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('CN'), 1), AtomID(pose.residue(1).atom_index('N'), 1), AtomID(pose.residue(1).atom_index('CA'), 1), AtomID(pose.residue(1).atom_index('C'), 1), i*0.0174532925)
		sc = sfxn3(pose)
		if sc < min_e:
			min_e = sc
			best = i

	pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('CN'), 1), AtomID(pose.residue(1).atom_index('N'), 1), AtomID(pose.residue(1).atom_index('CA'), 1), AtomID(pose.residue(1).atom_index('C'), 1), best*0.0174532925)


	min_e = 10**10; best = 0.0
	for i in range(0, 360):
		pose.conformation().set_torsion_angle(AtomID(pose.residue(res_length).atom_index('O'), res_length), AtomID(pose.residue(res_length).atom_index('C'), res_length), AtomID(pose.residue(res_length).atom_index('CA'), res_length), AtomID(pose.residue(res_length).atom_index('N'), res_length), i*0.0174532925)
		sc = sfxn3(pose)
		if sc < min_e:
			min_e = sc
			best = i

	pose.conformation().set_torsion_angle(AtomID(pose.residue(res_length).atom_index('O'), res_length), AtomID(pose.residue(res_length).atom_index('C'), res_length), AtomID(pose.residue(res_length).atom_index('CA'), res_length), AtomID(pose.residue(res_length).atom_index('N'), res_length), best*0.0174532925)
	pose.conformation().set_bond_angle( AtomID( pose.residue(1).atom_index('CN') , 1 ) , AtomID( pose.residue(1).atom_index('N') , 1 ) , AtomID( pose.residue(1).atom_index('CA') , 1 ), 119.122*0.0174533 )

	xsc = sfxn3(pose)
	for i in range(1, res_length+1):
		if pose.residue(i).has_property('N_METHYLATED'):
			min_e = 10**10; best = 0.0
			for j in range(0, 120):
				pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HN'), i), AtomID(pose.residue(i).atom_index('CN'), i), AtomID(pose.residue(i).atom_index('N'), i), AtomID(pose.residue(i).atom_index('CA'), i), j*0.0174532925)
				sc = sfxn3(pose)
				if sc < min_e:
					min_e = sc
					best = j

			pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HN'), i), AtomID(pose.residue(i).atom_index('CN'), i), AtomID(pose.residue(i).atom_index('N'), i), AtomID(pose.residue(i).atom_index('CA'), i), best*0.0174532925)
		if pose.residue(i).has_property('L_AA'):
			min_e = 10**10; best = 0.0
			for j in range(0, 120):
				pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HB'), i), AtomID(pose.residue(i).atom_index('CB'), i), AtomID(pose.residue(i).atom_index('CA'), i), AtomID(pose.residue(i).atom_index('N'), i), j*0.0174532925)
				sc = sfxn(pose)
				if sc < min_e:
					min_e = sc
					best = j

			pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HB'), i), AtomID(pose.residue(i).atom_index('CB'), i), AtomID(pose.residue(i).atom_index('CA'), i), AtomID(pose.residue(i).atom_index('N'), i), best*0.0174532925)
		
		if pose.residue(i).has_property('BETA_AA'):
			min_e = 10**10; best = 0.0
			for j in range(0, 120):
				pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HB'), i), AtomID(pose.residue(i).atom_index('CB'), i), AtomID(pose.residue(i).atom_index('CA'), i), AtomID(pose.residue(i).atom_index('N'), i), j*0.0174532925)
				sc = sfxn(pose)
				if sc < min_e:
					min_e = sc
					best = j

			pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HB'), i), AtomID(pose.residue(i).atom_index('CB'), i), AtomID(pose.residue(i).atom_index('CA'), i), AtomID(pose.residue(i).atom_index('N'), i), best*0.0174532925)


	#print 'before N-CH3 minimization', xsc
	#print 'after N-CH3 minimization', sfxn3(pose)
	fin_score = sfxn3(pose)
	if sys.argv[2] == 'SAR':
		path = '/suppscr/baker/twcraven/cyclic_peptides/' + 'poly_SAR_' + str(res_length) + 'mers/' + 'poly_SAR_' + str(res_length) + 'mer_' + str(fin_score) + '.pdb'
	elif sys.argv[2] == 'beta':
		path = '/suppscr/baker/twcraven/cyclic_peptides/' + 'poly_SAR_[B3A_res2]_' + str(res_length) + 'mers/' + 'poly_SAR_[B3A_res2]' + str(res_length) + 'mer_' + str(fin_score) + '.pdb'
	elif sys.argv[2] == 'alpha':
		path = '/suppscr/baker/twcraven/cyclic_peptides/' + 'poly_SAR_[ALA_res2]_' + str(res_length) + 'mers/' + 'poly_SAR_[ALA_res2]' + str(res_length) + 'mer_' + str(fin_score) + '.pdb'
	elif sys.argv[2] == 'oligourea':
		path = '/suppscr/baker/twcraven/cyclic_peptides/' + 'poly_SAR_[BB0_res2]_' + str(res_length) + 'mers/' + 'poly_SAR_[BB0_res2]' + str(res_length) + 'mer_' + str(fin_score) + '.pdb'

	N_coords = pose.residue(1).xyz(pose.residue(1).atom_index('N'))
	C_coords = pose.residue(len(pose.sequence())).xyz(pose.residue(len(pose.sequence())).atom_index('C'))
	
	if fin_score < 10.0 and dist(N_coords, C_coords) < 2.0:
		pose.dump_pdb(path)





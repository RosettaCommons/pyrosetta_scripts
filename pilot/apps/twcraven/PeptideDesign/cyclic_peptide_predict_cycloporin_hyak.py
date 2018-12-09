from pyrosetta import *
import numpy as np
import random
from rosetta.protocols.cyclic_peptide import DeclareBond
from rosetta.core.scoring.func import HarmonicFunc
import math

init('-database /suppscr/baker/twcraven/database -override_rsd_type_limit -ex1 -ex2 -ignore_unrecognized_res -mute all -output_virtual -bbdep_omega false -corrections::beta_nov16')

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
#help(sfxn)

movemap = MoveMap()
movemap.set_chi(True)
movemap.set_bb(True)
minmover = rosetta.protocols.simple_moves.MinMover(movemap, sfxn, "linmin", 0.001, True )

n_struct = int(sys.argv[1])

for k in range(0, n_struct):
	pose = Pose()
	res_length = 11
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY_p:N_Methylation_noHs')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation_noHs')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation_noHs')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('DALA')), True)
	
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation_noHs')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation_noHs')), True)
	
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation_noHs')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation_noHs')), True)
	pose.append_residue_by_bond(rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA')), True)
		
	"""	
	tx = []
	ty = []
	for i in range(0, 10000):
		torsions = pyrosetta.rosetta.utility.vector1_double()
		rama.random_mainchain_torsions(pose.conformation(), pose.residue(2).type(), pose.residue(3).type(), torsions)
		tx.append(torsions[1])
		ty.append(torsions[2])
		
	fig = plt.figure()
	ax = fig.add_subplot(111)

	#ax.scatter(tx, ty, '.', color='black', linewidth = 1)
	plt.scatter(tx, ty)

	plt.xlim(-180, 180)
	plt.ylim(-180, 180)
	plt.grid(True)
	plt.show()
	sys.exit(1)
	"""


	for i in range(1, res_length):
		torsions = pyrosetta.rosetta.utility.vector1_double()
		rama.random_mainchain_torsions(pose.conformation(), pose.residue(i).type(), pose.residue(i+1).type(), torsions)
		#print torsions
		

		if pose.residue(i+1).has_property('N_METHYLATED'):
			ri = np.random.randint(2, size=1)
			ri = 1
			if ri == 1:
				pose.set_omega(i, 180.0)
			else:
				pose.set_omega(i, 0.0)
		else:
			pose.set_omega(i, 180.0)
		
		for j in range(1, len(pose.residue(i).mainchain_torsions())):
			pose.set_torsion( rosetta.core.id.TorsionID(i, rosetta.core.id.BB, j), torsions[j] )
			
	pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('CN'), 1), AtomID(pose.residue(1).atom_index('N'), 1), AtomID(pose.residue(1).atom_index('CA'), 1), AtomID(pose.residue(1).atom_index('C'), 1), random.uniform(0, 360)*0.0174532925)
	pose.conformation().set_torsion_angle(AtomID(pose.residue(res_length).atom_index('O'), res_length), AtomID(pose.residue(res_length).atom_index('C'), res_length), AtomID(pose.residue(res_length).atom_index('CA'), res_length), AtomID(pose.residue(res_length).atom_index('N'), res_length), random.uniform(-140, 10)*0.0174532925)

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
	#print anchor_res, middle_loop_res
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
	#genkic.add_perturber('sample_cis_peptide_bond')
	#genkic.add_value_to_perturber_value_list( 2, 0.1)

	i = anchor_res + 1
	while (i!=anchor_res):
		if i > cyclization_point_start:
			i = cyclization_point_end
		elif i == anchor_res:
			i = i#do nothing
		else:
			genkic.add_residue_to_perturber_residue_list(1, i)
			#if i != res_length:
			#	if pose.residue(i+1).has_property('N_METHYLATED'):
			#		genkic.add_residue_to_perturber_residue_list(2, i)
		i += 1

	genkic.add_filter('loop_bump_check')

	ri = np.random.randint(2, size=1)
	ri = 1
	if ri == 1:
		genkic.close_bond(res_length, 'C', 1, 'N', 0, '', 0, '', 1.328, 116.206, 121.702, 180, False, False)

	else:
		genkic.close_bond(res_length, 'C', 1, 'N', 0, '', 0, '', 1.328, 116.206, 121.702, 0, False, False)

	genkic.set_pivot_atoms( first_loop_res, 'CA', middle_loop_res, 'CA', last_loop_res, 'CA' );
	genkic.set_closure_attempts(1000)
	genkic.set_selector_scorefunction(sfxn2)
	genkic.set_selector_type('lowest_energy_selector')
	genkic.apply(pose)

	#pose.replace_residue(6, rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY_p:N_Methylation')), True)
	pose.replace_residue(1,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('GLY_p:N_Methylation')), True)
	pose.replace_residue(2,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation')), True)
	pose.replace_residue(4,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation')), True)
	pose.replace_residue(7,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation')), True)
	pose.replace_residue(8,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation')), True)
	pose.replace_residue(9,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation')), True)
	pose.replace_residue(10,rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA_p:N_Methylation')), True)


	dhf = HarmonicFunc( 180.0*0.0174533 , 0.1 )
	dhc = rosetta.core.scoring.constraints.DihedralConstraint( AtomID(pose.residue(1).atom_index('N'), 1) , AtomID(pose.residue(res_length).atom_index('CA'), res_length) , AtomID(pose.residue(res_length).atom_index('C'), res_length) , AtomID(pose.residue(res_length).atom_index('O'), res_length) , dhf )
	pose.add_constraint( dhc )

	dhf = HarmonicFunc( 180.0*0.0174533 , 0.1 )
	dhc = rosetta.core.scoring.constraints.DihedralConstraint( AtomID(pose.residue(1).atom_index('CA'), 1) , AtomID(pose.residue(res_length).atom_index('C'), res_length) , AtomID(pose.residue(1).atom_index('N'), 1) , AtomID(pose.residue(1).atom_index('CN'), 1) , dhf )
	pose.add_constraint( dhc )


	min_e = 10**10; best = 0.0
	for i in range(0, 360):
		pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('CN'), 1), AtomID(pose.residue(1).atom_index('N'), 1), AtomID(pose.residue(1).atom_index('CA'), 1), AtomID(pose.residue(1).atom_index('C'), 1), i*0.0174532925)
		sc = sfxn(pose)
		if sc < min_e:
			min_e = sc
			best = i

	pose.conformation().set_torsion_angle(AtomID(pose.residue(1).atom_index('CN'), 1), AtomID(pose.residue(1).atom_index('N'), 1), AtomID(pose.residue(1).atom_index('CA'), 1), AtomID(pose.residue(1).atom_index('C'), 1), best*0.0174532925)


	min_e = 10**10; best = 0.0
	for i in range(0, 360):
		pose.conformation().set_torsion_angle(AtomID(pose.residue(res_length).atom_index('O'), res_length), AtomID(pose.residue(res_length).atom_index('C'), res_length), AtomID(pose.residue(res_length).atom_index('CA'), res_length), AtomID(pose.residue(res_length).atom_index('N'), res_length), i*0.0174532925)
		sc = sfxn(pose)
		if sc < min_e:
			min_e = sc
			best = i

	pose.conformation().set_torsion_angle(AtomID(pose.residue(res_length).atom_index('O'), res_length), AtomID(pose.residue(res_length).atom_index('C'), res_length), AtomID(pose.residue(res_length).atom_index('CA'), res_length), AtomID(pose.residue(res_length).atom_index('N'), res_length), best*0.0174532925)
	pose.conformation().set_bond_angle( AtomID( pose.residue(1).atom_index('CN') , 1 ) , AtomID( pose.residue(1).atom_index('N') , 1 ) , AtomID( pose.residue(1).atom_index('CA') , 1 ), 119.122*0.0174533 )

	xsc = sfxn(pose)
	for i in range(1, res_length+1):
		if pose.residue(i).has_property('N_METHYLATED'):
			min_e = 10**10; best = 0.0
			for j in range(0, 120):
				pose.conformation().set_torsion_angle(AtomID(pose.residue(i).atom_index('1HN'), i), AtomID(pose.residue(i).atom_index('CN'), i), AtomID(pose.residue(i).atom_index('N'), i), AtomID(pose.residue(i).atom_index('CA'), i), j*0.0174532925)
				sc = sfxn(pose)
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

	N_coords = pose.residue(1).xyz(pose.residue(1).atom_index('N'))
	C_coords = pose.residue(len(pose.sequence())).xyz(pose.residue(len(pose.sequence())).atom_index('C'))

	fin_score = sfxn(pose)
	#print fin_score
	if fin_score < 75.0 and dist(N_coords, C_coords) < 2.0:
		pose.dump_pdb('/suppscr/baker/twcraven/cyclosporinA_backbones/cyclosporinA_'+str(sfxn(pose))+'.pdb')



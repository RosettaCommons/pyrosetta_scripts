#!/home/sheffler/venv/david/bin/python
from sys import argv
from pyrosetta import *
from pdb_utils_noclass import *
#from pyrosetta.toolbox import pose_from_rcsb
from repeat_utils import *
import random,string
from math import *
from xyzMath import *
#from two_sided_design_pyr import *
import pickle
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from pyrosetta.rosetta.numeric import xyzVector_double_t as V3
from pyrosetta.rosetta.numeric import xyzMatrix_double_t as M3
from hash_subclass import *
import argparse
import numpy as np
import itertools
import os

def get_argparse():
    parser = argparse.ArgumentParser(description='Dock de novo generated peptide backbones against repeat proteins')
    parser.add_argument('input_file',
                   help='list of repeat proteins')
    parser.add_argument('--dock_res', type=int, dest='nbins',
                   default=6,
                   help='peptide-repeat protein docking resolution (higher number is higher res and slower)')
    parser.add_argument('--npeptides', type=int, dest='npept',
                   default=200,
                   help='number of peptides to generate')
    parser.add_argument('--phi_range', type=float, dest='phi_range', nargs=2,
                   default=[-180.,180.],
                   help='lower and upper limits of peptide phi distribution')
    parser.add_argument('--psi_range', type=float, dest='psi_range', nargs=2,
                   default=[-180.,180.],
                   help='lower and upper limits of peptide psi distribution')
    parser.add_argument('--phi_step', type=float, dest='phi_step', nargs=1,
                   default=10,
                   help='')
    parser.add_argument('--psi_step', type=float, dest='psi_step', nargs=1,
                   default=10,
                   help='')
    parser.add_argument('--Nrepeats', type=int, dest='Nrepeats', 
                   default=6,
                   help='number of repeat units in peptide')
    parser.add_argument('--repeat_length', type=int, dest='repeat_length', 
                   default=2,
                   help='length of repeat unit')
    parser.add_argument('--skip_docking',type=bool, dest='skip_docking',
                   default=False,                    
                   help='skip docking step')
    parser.add_argument('--use_hash', type=bool, dest='use_hash',
                   default=True,
                   help='use bidentate hbond hash (or other geometric property) to filter docks ')
    parser.add_argument('--hash_file_name', dest='hash_file', 
                   default='bb_asn_dict_combo_1.0_rot0',
                   help='name of hash file (required if --use_hash is specified if different from default) ')
    parser.add_argument('--repeat_match_trans', type= float, dest='trans_threshold', 
                   default='2.0',
                   help='maximum difference in helical translation per repeat (Angstroms) between peptide and repeat protein')
    parser.add_argument('--repeat_match_rot', type=float, dest='rot_threshold', 
                   default='1.0',
                   help='maximum difference in helical rotation per repeat (radians) between peptide and repeat protein')
    parser.add_argument('--use_lineal_scan', type=bool, dest='use_lineal_scan',
                   default=False,
                   help='scan all the angles combinations')
    return parser

def choose_torsions_randomly(phi_range,
                             psi_range,
                             aa_name,
                             rama_score_limit=1.0,
                             MAX_TRIES=20):
    #MAX_TRIES=20
    #ntries=0
    #found=0
    sfx_rama = rosetta.core.scoring.Ramachandran()
    res_type = rosetta.core.chemical.aa_from_name(aa_name)
    torsions_list=[]
    for i in xrange(MAX_TRIES):
      #ntries=ntries+1
      phi = random.uniform(phi_range[0],phi_range[1])
      psi = random.uniform(psi_range[0],psi_range[1])
      if sfx_rama.eval_rama_score_residue(res_type,phi,psi) < 1.0:
        #found=1
	print "Rand torsions found: %0.3f, %0.3f"%(phi,psi)
        torsions_list.append([phi,psi])

    if len(torsions_list) > 0:
        return torsions_list
    else:
        return False
    #return [phi, psi]  # if not found, then return last phi psi pair sampled


def choose_torsions_from_linear_scan(phi_range=[-180.,180.],
                                     psi_range=[-180.,180.],
                                     phi_step=5.0,
				     psi_step=5.0,
                                     aa_name="ASN",
                                     rama_score_limit=0.0
                                     ):
    sfx_rama = rosetta.core.scoring.Ramachandran()
    res_type = rosetta.core.chemical.aa_from_name(aa_name)
    torsions_list=[]
    for phi_i in np.arange(phi_range[0],phi_range[1]+phi_step,phi_step):
        for psi_j in np.arange(psi_range[0],psi_range[1]+phi_step,psi_step):
            if sfx_rama.eval_rama_score_residue(res_type,phi_i,psi_j) < rama_score_limit:
                torsions_list.append([phi_i,psi_j])
    if len(torsions_list) > 0:
        return torsions_list
    else:
        return False

def switch_pose_to_centroid(in_pose):
  switchToCentroid = SwitchResidueTypeSetMover("centroid")
  switchToCentroid.apply(in_pose)
  return(in_pose)
 
def generate_peptides(repeat_L,
                      N_repeats,
		      N_peptides,
                      phi_range,
                      psi_range,
                      phi_step=5.0,
                      psi_step=5.0,
                      aa_for_rama_space="ASN",
                      rama_score_limit=0.0,
                      fakeResForBuildingPeptide="A",
                      centroid_VDWscore_perRes_limit=0.0,
                      use_lineal_scan=False,
                      dump_pdb=False):
  peptide_len=N_repeats*repeat_L
  print "Generating parameters for a %d a.a. peptide with NumRep=%d, LenRep=%d"%(peptide_len, N_repeats, repeat_L)
  torsions=[]
  if use_lineal_scan:
      print "Using lineal torsion scan"
      torsions=choose_torsions_from_linear_scan(phi_range=phi_range,
                                                psi_range=psi_range,
                                                phi_step=phi_step,
                                                psi_step=psi_step,
                                                aa_name=aa_for_rama_space,
                                                rama_score_limit=rama_score_limit)
  else:
      print "Using random torsion scan"
      torsions=choose_torsions_randomly(phi_range=phi_range,
                                        psi_range=psi_range,
                                        aa_name=aa_for_rama_space,
                                        rama_score_limit=rama_score_limit,
                                        MAX_TRIES=20)

  print "Number of found allowed torsions:", len(torsions)

  p=Pose()
  pdbs={}
  peptide_len=N_repeats*repeat_L
  pyrosetta.rosetta.core.pose.make_pose_from_sequence(p, fakeResForBuildingPeptide*peptide_len,"fa_standard")

  #Wiggling the pose around is faster in centroid mode
  switch_pose_to_centroid(p)

  for ires in xrange(peptide_len):
       p.set_omega(ires+1,180.0)

  scorefxn_cen = pyrosetta.create_score_function("cen_std") #pyrosetta.get_fa_scorefxn() #create_score_function("cen_std") 
  print scorefxn_cen

  helical_params=[]
  phi_psiL=[]

  print "Target number of phi_psi combinations: %d"% len(list(itertools.product(torsions, repeat=repeat_L)))
  count=0
  random.shuffle(torsions) #Randomnize this in case we don't need all of them
  print torsions
  for phi_psi_product in itertools.product(torsions, repeat=repeat_L):
      count+=1;
      if (count%1e4==0):
	print "Working on combination:",  count
      for i_pair in xrange(repeat_L):
         for j_repeat in xrange(N_repeats):
             #print (i_pair+(j_repeat*repeat_L))
             p.set_phi((i_pair+(j_repeat*repeat_L))+1,phi_psi_product[i_pair][0])
             p.set_psi((i_pair+(j_repeat*repeat_L))+1,phi_psi_product[i_pair][1])
      #break
      scorefxn_cen(p)
      rep_score=p.energies().total_energies().get(pyrosetta.rosetta.core.scoring.ScoreType.vdw)
      if (rep_score/peptide_len) <= centroid_VDWscore_perRes_limit:
             #Use residues 1,3 to determine the peptide parameters, meaning that the min peptide size is 3a.a. 
             #Note:(this might not be the right thing to do, this about helixes that wrap around)
             res1=[Vec(p.residue(1).xyz("N")),Vec(p.residue(1).xyz("CA")),Vec(p.residue(1).xyz("C"))]
             res2=[Vec(p.residue(3).xyz("N")),Vec(p.residue(3).xyz("CA")),Vec(p.residue(3).xyz("C"))]

             #It seems get_helix_params might fail in certain scenarios, for example:
             #res1=[Vec( 2.228983, -0.683234, -3.163434 ), Vec( 2.582034, 0.000000, -1.924760 ), Vec( 1.462372, 0.926429, -1.468252 )]
             #res2=[Vec( -1.607544, 2.461074, 0.493442 ), Vec( -1.958726, 2.462015, 1.908517 ), Vec( -0.837940, 1.871868, 2.754614 )]
             try: #Will use this try to skip the bad ones (like the one in the example before)
                 trans, radius, ang = get_helix_params( res1, res2 )
                 #helical_params.append( (trans,radius,ang) )
                 p=center_on_z_axis(res1,res2,p)
                 #pdbs[count]=p.clone()
                 #if dump_pdb:   
                 #    p.dump_pdb('testpdbs/test_%05d_tf.pdb'%count)
                 #print count, rep_score, phi_psi_product, trans, radius, ang
                 yield phi_psi_product, p.clone(), (trans,radius,ang)
             except:
                 print "That was an incommesurable peptide!!! Will skipp it!"
                 continue
             #yield phi_psi_product, p.clone(), (trans,radius,ang)

def dock_peptide(p,nbins): 
# have match between repeat protein q and peptide p. now sample helical degrees of freedom to dock p on q
        for dis in range(nbins):
         for angle in range(nbins*6):
          deg=360.*float(angle)/float(nbins*6.)
          dist=-3.+6*(float(dis)/float(nbins))
#don't want p to be moved. the following is a bit silly          
          
          b=rotate_around_z( Vec(0.,0.,1.0),deg,p.clone())
          c=translate_along_z(dist,b.clone())
#          c.append_pose_by_jump(q.clone(),1)
#          c.dump_pdb('%s_%s_combo_%s_%s_tf.pdb'%(DHR[0:base],p_index,angle,dist))
          yield(c, deg,dist)
  
def compare_params(pept_gen,
                   DHR,
                   names, 
                   trans_threshold, 
                   rot_threshold,
                   max_matches,
                   b_dump_pdb=False):
    matches={}
    for d in names:
        matches[d]=[]
    #for j in range(len(pept)):

    Nmatch=0
    for torsions, p, params in pept_gen: #yield phi_psi_product, p.clone(), (trans,radius,ang)
        for i in range(len(DHR)):
            d=DHR[i]
            for n in range(1,4): # check for 1 - 2 peptide repeats equivalent to 1 repeat protein repeat
             if ( (abs(d[0] -n*params[0]) < trans_threshold) and (abs(d[2]-n*params[2]) < rot_threshold)) :
                matches[names[i]].append( (p.clone(),params,torsions) )
		print "Found one match (#%d):"%len(matches[names[i]]), matches[names[i]][-1][1], matches[names[i]][-1][2]
		if b_dump_pdb: 
                    p.dump_pdb('testpdbs/test_%05d_match.pdb'%Nmatch)
                Nmatch+=1
                ##break #Debug
	if (Nmatch>=max_matches):
            break
    return Nmatch,matches

def main():
   ############################################################
   parser = get_argparse()
   args=parser.parse_args()
   init_pyrosetta()
   
   
   #input_file=argv[1] # list of repeat proteins
   ##nbins=int(argv[2]) # resolution of docking grid sampling
   #npept=int(argv[3])
   #Nstruct, angle variance, output_pdb

   print 'get repeat protein params '
   DHR_params, DHR_arcs, names, lengths, rep_structs = calc_repeat_protein_params_ws(args.input_file)
   
   # use generator to avoid memory cost of storing all structures
   print 'set up  peptide backbone generator '
   print args.repeat_length, args.Nrepeats, args.npept, args.phi_range
   pept_gen = generate_peptides(args.repeat_length,
                                args.Nrepeats,
                                args.npept,
                                phi_range=args.phi_range,
                                psi_range=args.psi_range,
				phi_step=args.phi_step,
                                psi_step=args.psi_step,
                                aa_for_rama_space="ASN",
                                rama_score_limit=0.0,
                                fakeResForBuildingPeptide="A",
                                centroid_VDWscore_perRes_limit=0.0,
                                use_lineal_scan=(args.use_lineal_scan==True),
                                dump_pdb="False"
                                )  #repeat_length, Nstruct, angle variance, bool_output_pdb
   
   print 'generate peptides and compare helical params to those of repeat proteins'
   Nmatch, matches=compare_params(pept_gen,
                                  DHR_params,
                                  names,
                                  args.trans_threshold,
                                  args.rot_threshold,
                                  max_matches=args.npept,
                                  b_dump_pdb=True ) #Return: p.clone(),params,torsions
   print('Number of matches: %s '%Nmatch)
   
   for DHR in matches.keys():
      #print DHR, DHR_arcs[DHR],len(matches[DHR])
      q=rep_structs[DHR]
      pept_match=[]
      for i in range(len(matches[DHR])):
         print 'docking peptide %s'%i
         p=matches[DHR][i][0]

         pept_params=matches[DHR][i][1]
         torsions=matches[DHR][i][2]
         arc_length=sqrt( pept_params[0]**2 + (pept_params[1] * sin(pept_params[2] )**2  ) )
         print('%s %.2f %.2f %.2f  %.2f'%(i,pept_params[0],pept_params[1],pept_params[2],arc_length))
         
         if (args.skip_docking==True): 
             print args.skip_docking
             continue
         dock_gen=dock_peptide(p,args.nbins)
         print 'evaluate number of contacts and clashes for each dock of peptide'
         good_matches,contact_hist, ntries=eval_contacts_pair(q, dock_gen, 20)  #need to fix now that have generator
         
         
            
         print 'Trial %d, matches: %d '%(ntries, len(good_matches))
         for match in good_matches:
            new_pdb=match[0].clone()
            n_sc_bb=0
            if (args.use_hash==True):
                  hash_nc=use_hash_rot(1.0,3.0,args.hash_file)
                  n_sc_bb=hash_nc.count_asn_bb(new_pdb,q)
                  print('number asn-bb hbonds %s'%n_sc_bb)
            if (args.use_hash==False) or n_sc_bb > 0:
               new_pdb=pdb.clone() 
               new_pdb.append_pose_by_jump(q.clone(),1)
               trans=round(pept_params[0],1)
               twist=round(pept_params[2],1)
               angle=round(match[1],2)
               dist=round(match[2],2)

               out_name="testpdbs/%s_%s_%s_%0.2f_%0.2f_%0.2f.pdb"%(os.path.basename(DHR)[:-4],
                                                                n_sc_bb,
                                                                str(torsions).replace(" ", "").replace("[", "").replace("]", "").replace(",", ""),
                                                                angle,
                                                                dist,
                                                                twist)
               print "OUT MaTcH:", out_name
               new_pdb.dump_pdb(out_name)

   

if __name__ == "__main__":
    main()

#!/software/miniconda3/envs/pyrosetta3/bin//python
from sys import argv
from pyrosetta import *
from pdb_utils_noclass import *
from pyrosetta.toolbox import pose_from_rcsb
from repeat_utils_w_contact_list import *
import random,string
from math import *
from rif.legacy.xyzMath import *
#from two_sided_design_pyr import *
import pickle
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from rosetta.numeric import xyzVector_double_t as V3
from rosetta.numeric import xyzMatrix_double_t as M3
from hash_subclass import *
import argparse

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
                   default=[-180.,189.],
                   help='lower and upper limits of peptide psi distribution')
    parser.add_argument('--Nrepeats', type=int, dest='Nrepeats',
                   default=6,
                   help='number of repeat units in peptide')
    parser.add_argument('--repeat_length', type=int, dest='repeat_length',
                   default=2,
                   help='length of repeat unit')
    parser.add_argument('--skip_docking',type=int, dest='skip_dock',
                   default='0',
                   help='skip docking step')
    parser.add_argument('--use_hash', type=int, dest='use_hash',
                   default='1',
                   help='use bidentate hbond hash (or other geometric property) to filter docks ')
    parser.add_argument('--hash_file_name', dest='hash_file',
                   default='bb_asn_dict_combo_1.0_rot0',
                   help='name of hash file (required if --use_hash is specified if different from default) ')
    parser.add_argument('--repeat_match_trans', type= float, dest='trans_threshold',
                   default='1.5',
                   help='maximum difference in helical translation per repeat (Angstroms) between peptide and repeat protein')
    parser.add_argument('--repeat_match_rot', type=float, dest='rot_threshold',
                   default='0.2',
                   help='maximum difference in helical rotation per repeat (radians) between peptide and repeat protein')
    parser.add_argument('--struct_dir', type=str, dest='struct_dir',
                   default="/work/baker/repeat_peptide/designs/",
                   help='This is where the repeat protein designs/xtals pdbs are located')
    return parser

def choose_torsions(phi_range,psi_range,aa_type):
    MAX_TRIES=20
    ntries=0
    found=0
    for i in range(MAX_TRIES):
      ntries=ntries+1
      phi = random.uniform(phi_range[0],phi_range[1])
      psi = random.uniform(psi_range[0],psi_range[1])
      if eval_rama(aa_type,phi,psi) < 1.0:
        found=1
        break

    if found==0: print('no torsions found')
    return phi,psi  # if not found, then return last phi psi pair sampled

def generate_peptides(repeat_L,N_repeats,Nsamples,phi_range,psi_range,dump_pdb):
  p=Pose()
  pdbs={}
  length=N_repeats*repeat_L
  seq=""
  for i in range(length):seq=seq+"A"
  rosetta.core.pose.make_pose_from_sequence(p, seq,"fa_standard")
  torsions=[]
#  helical_params=[]
  for i in range(Nsamples):
    torsions=[]
    for j in range(repeat_L):
        phi,psi=  choose_torsions(phi_range,psi_range,'THR')
        torsions.append( (phi,psi) )


    for start in [1+repeat_L*n for n in range(N_repeats)]:

      for k in range(repeat_L):
        p.set_phi(start+k,torsions[k][0])
        p.set_psi(start+k,torsions[k][1])
        p.set_omega(start+k,180.)

    tor_list=[]
    for torsion in torsions:
        tor_list.append(torsion[0])
        tor_list.append(torsion[1])

    if dump_pdb:    p.dump_pdb('test_%s.pdb'%i)

    res1=[Vec(p.residue(1).xyz("N")),Vec(p.residue(1).xyz("CA")),Vec(p.residue(1).xyz("C"))]
    res2=[Vec(p.residue(repeat_L+1).xyz("N")),Vec(p.residue(repeat_L+1).xyz("CA")),Vec(p.residue(repeat_L+1).xyz("C"))]
    trans, radius, ang =  get_helix_params(res1,res2)
#    helical_params.append( (trans,radius,ang) )
    p=center_on_z_axis(res1,res2,p)
    pdbs[i]=p.clone()
    if dump_pdb:   p.dump_pdb('test_%s_tf.pdb'%i)
    yield tor_list, p.clone(), (trans,radius,ang)
  # return torsions, pdbs, helical_params

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

def compare_params(pept_gen, DHR,names, trans_threshold, rot_threshold):
#    trans_threshold=2.0  # difference in translation/repeat (Angstroms)
#    rot_threshold=0.6  # difference in rotation/repeat (radians)
    matches={}
    Nmatch=0
    for d in names:
        matches[d]=[]
    #for j in range(len(pept)):
    for torsions, pdb, p in pept_gen:
	#print pdb,pdb.phi(4),pdb.psi(4)
	#print torsions
        for i in range(len(DHR)):
            d=DHR[i]
            for n in range(1,4): # check for 1 - 2 peptide repeats equivalent to 1 repeat protein repeat
             if ( abs(d[0] -n*p[0]) < trans_threshold and abs(d[2]-n*p[2]) < rot_threshold) :
                #print('match  %s %s'%(names[i],j))
                matches[names[i]].append( (pdb.clone(),p,torsions) )
                Nmatch=Nmatch+1
                break
    return Nmatch,matches


############################################################
parser = get_argparse()
args=parser.parse_args()
init_pyrosetta()
#input_file=argv[1] # list of repeat proteins
##nbins=int(argv[2]) # resolution of docking grid sampling
#npept=int(argv[3])
#Nstruct, angle variance, output_pdb

# use generator to avoid memory cost of storing all structures
print('set up  peptide backbone generator ')
pept_gen = generate_peptides(args.repeat_length,args.Nrepeats,args.npept,args.phi_range,args.psi_range,0)  #repeat_length, Nstruct, angle variance, output_pdb

print('get repeat protein params ')
DHR_params, DHR_arcs, names, lengths, rep_structs = calc_repeat_protein_params_ws(args.input_file,args.struct_dir,offset=0)
#DHR_params, DHR_arcs, names, lengths, rep_structs = calc_repeat_protein_params_ws(args.input_file,args.struct_dir,offset=10)
print('generate peptides and compare helical params to those of repeat proteins')
Nmatch, matches=compare_params(pept_gen,DHR_params,names,args.trans_threshold,args.rot_threshold)
print(('Number of matches: %s '%Nmatch))

if args.use_hash: hash_nc=use_hash_rot(1.0,3.0,args.hash_file)
for DHR in list(matches.keys()):
    print(DHR, DHR_arcs[DHR],len(matches[DHR]))
    q=rep_structs[DHR]
    pept_match=[]
    for i in range(len(matches[DHR])):
        print('docking peptide %s'%i)
        p=matches[DHR][i][0]
        pept_params=matches[DHR][i][1]
        torsions=matches[DHR][i][2]
        arc_length=sqrt( pept_params[0]**2 + (pept_params[1] * sin(pept_params[2] )**2  ) )
        str='%s %.2f %.2f %.2f  %.2f '%(i,pept_params[0],pept_params[1],pept_params[2],arc_length)
        for tor in torsions: str+=' %.2f '%tor
        print(str)
#        print('%s %.2f %.2f %.2f  %.2f '%(i,pept_params[0],pept_params[1],pept_params[2],arc_length))
#        print(torsions)
        if '.' in DHR:
            base=DHR.index('.')
        else:
            base=len(DHR)

        if (args.skip_dock): continue
        dock_gen=dock_peptide(p,args.nbins)
        print('evaluate number of contacts and clashes for each dock of peptide')
        good_matches,contact_hist, ntries=eval_contacts_pair(q,dock_gen,20)  #need to fix now that have generator

        print('contacts: ', ntries, len(good_matches),contact_hist)
        for match in good_matches:
         pdb=match[0].clone()
         contact_set=match[4]
         n_sc_bb=0
         if args.use_hash:
#  screen matches to count number of asn-bb hbonds
             n_sc_bb=hash_nc.count_asn_bb_w_contact_set(pdb,q,contact_set)
             print(('number asn-bb hbonds %s'%n_sc_bb))
         if not(args.use_hash) or n_sc_bb > 0:
            new_pdb=pdb.clone()
            new_pdb.append_pose_by_jump(q.clone(),1)
 #           params=match[1]
#            p_index=params[0]
            trans=round(pept_params[0],1)
            twist=round(pept_params[2],1)
            angle=round(match[1],2)
            dist=round(match[2],2)

            torsions_string ='_'.join(['%4.1f' % (t,) for t in torsions])

            new_pdb.dump_pdb('%s_%s_%s_%s_%s.pdb'%(DHR[0:base],n_sc_bb,torsions_string,angle,dist))
            print(('%s_%s_%s_%s_%s.pdb'%(DHR[0:base],n_sc_bb,torsions_string,angle,dist)))

            torsions_string_space =' '.join(['%4.1f' % (t,) for t in torsions])
            print(('%s %s %s %s %s %s %s %s'%(DHR[0:base],angle,dist,torsions_string_space,match[2],n_sc_bb,trans,twist)))

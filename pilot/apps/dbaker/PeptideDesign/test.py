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
    parser.add_argument('--psi_range', type=float, dest='phi_range', nargs=2,
                   default=[-180.,189.],
                   help='lower and upper limits of peptide psi distribution')
    parser.add_argument('--Nrepeats', type=int, dest='Nrepeat', 
                   default=6,
                   help='number of repeat units in peptide')
    parser.add_argument('--repeat_length', type=int, dest='repeat_length', 
                   default=2,
                   help='length of repeat unit')
    parser.add_argument('--skip_docking',type=bool, dest='skip_dock',
                   default='0',                    
                   help='skip docking step')
    parser.add_argument('--use_hash', type=bool, dest='use_hash',
                   default='1',
                   help='use bidentate hbond hash (or other geometric property) to filter docks ')
    parser.add_argument('--hash_file_name', dest='hash_file', 
                   default='bb_asn_dict_combo_1.0_rot0',
                   help='name of hash file (required if --use_hash is specified) ')
    
    
 
    
    
    
    return parser

parser = get_argparse()
args=parser.parse_args()
print args
print args.input_file,args.nbins

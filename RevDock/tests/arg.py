
import argparser 

def get_options():

    parser = argparse.ArgumentParser(description='Version: 1.3.1')
    parser.add_argument('-R', metavar='MGL root', nargs='?', default=std_in,
                        help='MGL_ROOT environmental variable - path to mgltools')
    parser.add_argument('-l', metavar='ligand.pdbqt', nargs='?', default=std_in,
                        help='path to the ligand.pdbqt file')
    parser.add_argument('-d', metavar='receptor directory', required=True,
                        help='path to the reference genome fasta file')
    parser.add_argument('-A', metavar='annotation', required=True,
                        help='"grch37" (GENCODE V24lift37 canonical annotation file in '
                             'package), "grch38" (GENCODE V24 canonical annotation file in '
                             'package), or path to a similar custom gene annotation file')
    parser.add_argument('-D', metavar='distance', nargs='?', default=50,
                        type=int, choices=range(0, 5000),
                        help='maximum distance between the variant and gained/lost splice '
                             'site, defaults to 50')
    parser.add_argument('-M', metavar='mask', nargs='?', default=0,
                        type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and '
                             'unannotated acceptor/donor loss, defaults to 0')
    args = parser.parse_args()

    return args
[DOCKING]
### Configuration file for the reverse docking procedure ###
### change file names according with your paths ####

# MGL_ROOT environmental variable - path to mgltools   
MGL_ROOT = /home/marco/MGLTools-1.5.6

# Specific input/output files
ligand = ./Docking_data/HTML_PLP.pdbqt            # ligand file
receptor_dir = ./Mouse/pdb                        # path to receptor pdb files
coord_file = ./Mouse/Mouse_Orf_pred_mini.coord          # table of gridcenter:  receptor_name1 x y z
outdir = ./Mouse/outdir_HTML_PLP                  # directory of results (if exist, add progressive number)
processes=5                                       # number of processes of multiprocessing

[GPF]
npts=55,55,55 

[DPF]
ga_pop_size=200
ga_num_evals=1000000
ga_run=5

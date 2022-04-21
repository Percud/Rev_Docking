import re, sys, os, subprocess, glob, configparser, pandas as pd, numpy as np
from multiprocessing import Pool
from biopandas.pdb import PandasPdb
from Bio.PDB.vectors import *
import math
from pymol import cmd
cwd = os.getcwd()

#####    config parser    ####
config = configparser.ConfigParser(inline_comment_prefixes="#")
config.optionxform = str
try:
    config.read(sys.argv[1])
    locals().update(dict(config.items('DOCKING')))
    locals().update(dict(config.items('GPF')))
    locals().update(dict(config.items('DPF')))
except:
    print(
'''
rd_conf.txt is required as argument, please provide it.
USAGE : rev_docking.py rd_conf.txt
Example from RevDock/modules/tests/rd_conf.txt:
[DOCKING]
### Configuration file for the reverse docking procedure ###
###
### change file names according with your paths ####
###

# Specific 'input/output files'
all_site = False
ligand = /hpc/home/marco.malatesta1/Rev_Docking/Docking_data/HTML_PLP.pdbqt            # ligand file
receptor_dir = .                        # path to receptor pdb files
coord_file = ./Mouse_coord.csv        # table of gridcenter:  receptor_name1 x y z
processes = 60                                                  #                       receptor_name2 x y z
outdir = ./HTML_PLP_Mouse_flex                  # directory of results
[GPF]
npts=25 25 25  ## unit in Ångström
spacing=0.2    ## unit in Ångström

[DPF]
ga_num_evals=50000000  # number of energy evaluations
ga_run=200   # number of independent docking runs

'''
    )
    sys.exit(1)

####   make outdir   ####

i = 1
while os.path.exists(outdir+'_%s' % i):
    i += 1
out_dir = outdir+'_%s' % i+'/'
os.makedirs(out_dir)
if os.path.exists(outdir+'_%s' % i):
    print('writing output in '+out_dir)

####    read coordinates    ####

try:
    print(f'Reading {coord_file} ...')
    
    if eval(all_site):
        coordinates = (pd.read_csv(coord_file, sep='\t', comment='#').dropna().groupby('pdb')
                           .apply(lambda r: r[['x','y','z', 'lys', 'chain']].values.tolist())
                           .to_dict())
        print(f'Proceeding with all site docking')
        
    else:
        coordinates = (pd.read_csv(coord_file, sep='\t', comment='#').dropna().sort_values('chain').drop_duplicates('pdb').groupby('pdb')
                         .apply(lambda r: r[['x','y','z', 'res', 'chain']].values.tolist())
                         .to_dict())
        print(f'Proceeding with one site docking')
        
        
except:
    print(coord_file+" doesn't exist")
    sys.exit(1)

####   copy the ligand   ####

try:
    os.system(f'cp {ligand} {out_dir}')
    lig_name = os.path.splitext(os.path.basename(ligand))[0]
except:
    print(ligand+" doesn't exist")
    sys.exit(1)


    
####    reverse docking    ####

for pdb in coordinates:
    try:
        xyz = coordinates.get(pdb)
        pdb_name = os.path.splitext(pdb)[0]
        for c in xyz:
            i=1
            while os.path.exists(f'{out_dir}{pdb_name}_{i}.pdbqt'):
                i += 1 # incrementing file number
                
            if 'rec_files_dir' in locals():
                os.system(f'cp {rec_files_dir}/{pdb_name}_{i}.trg {out_dir}')
                os.system(f'cp {rec_files_dir}/{pdb_name}_{i}.pdbqt {out_dir}')
            else:

                ####    prepare receptor    ####
                print(re.sub('HETATM', '#HETATM', open(f'{receptor_dir}/{pdb}', 'r').read()), 
                      file=open(out_dir+pdb, 'w'))                   # comment HETATM in pdb file
                os.system(f"prepare_receptor -A 'checkhydrogens' -e -r {out_dir}{pdb} -o {out_dir}{pdb_name}_{i}.pdbqt") 

                ####    AGFR    ####
                os.chdir(out_dir)
                os.system(f'agfr -r {pdb_name}_{i}.pdbqt -b user {" ".join(map(lambda x: str(round(x,3)),c[0:3]))} {npts} -s {spacing} -f {c[4]}:LYS{str(int(c[3]))}')
                os.chdir(cwd)

            ####    ADFR    ####
            os.chdir(out_dir)
            os.system(f'adfr -t {pdb_name}_{i}.trg -l {os.path.basename(ligand)} -J {pdb_name}_{i} -e {ga_num_evals} -n {ga_run} -f -c {processes} -D 3')
            os.chdir(cwd)
        
    except:
        print(f'{pdb} not found')
        os.chdir(cwd)

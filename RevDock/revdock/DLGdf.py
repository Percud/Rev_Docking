####    make table from dlg directory with different kind of energy    #### 
####    USAGE: ./get_atom_energy.py MY_RES_IN_LIGAND MY_DLG_DIRECTORY

import glob, re, sys, os, math, pandas as pd, argparse
from Bio.PDB.vectors import *
from biopandas.pdb import PandasPdb

def distance(a,b):
    '''Return max distance given arrays in xyz of points'''
    x = [a[0], b[0]]
    y = [a[1], b[1]]
    z = [a[2], b[2]]
    return (((max(x)-min(x))**2+(max(y)-min(y))**2+(max(z)-min(z))**2)**0.5)

parser = argparse.ArgumentParser(prog = 'DLGdf',
                                 epilog = '''
                                            Returns a Dataframe.tsv given a three-letter code of the ligand and a dlg directory
                                            USAGE : DLGdf three_letter_code dlg_directory
                                            Example : DLGdf -i *.dlg -r GLU -d 2''')
parser.add_argument('-i', '--input',
                    nargs = '+',
                    required = True,
                    help = 'Define DLG input file')
parser.add_argument('-r', '--residue',
                    help = 'Define residue three-letter code to compute energy separately')
parser.add_argument('-d', '--distance',
                    help = 'Define atom number to measure distance from grid center')
parser.add_argument('-c', '--catalytic',
                    nargs = '+',
                    help = 'Define atom number to measure distance from grid center')
parser.add_argument('-a', '--angle',
                    nargs = '+',
                    help = 'Define atom number to measure dihdral')
args = parser.parse_args()

parser.add_argument('-o', '--output',
                    default = f"{os.path.dirname(args.input[0])}/{os.path.basename(os.path.dirname(args.input[0]))}.tsv",
                    help = 'Define output filename of tsv format table')
args = parser.parse_args()


dlg_dir = os.path.dirname(os.path.abspath(args.input[0]))  ## DLGs directory    

## print header of table

print('dlg file', 'Ligand', 
      'BCB (kcal/mol)', 'Run', 'LCB (kcal/mol)', 'Run', 
      'BCaaB (kcal/mol)', 'Run', 'LCaaB (kcal/mol)', 'Run', 
      'BCM (kcal/mol)', 'Runs', 'LCM (kcal/mol)', 'Runs', 
      'BCaaM (kcal/mol)', 'Runs', 'LCaaM (kcal/mol)', 'Runs', 
      'LC', 'Num in LC', '1LC/2LC %', 
      'Distance (A)', 'Dihedral',
      'Catalytic Residue', 'Catalytic Residue Occurrence', 'Catalytic LC',
      sep = '\t',
      file=open(f"{dlg_dir}/{os.path.basename(dlg_dir)}.tsv", 'w')) # tab headers

for dlg in args.input:
    try:
        tab, around, count = [], [], []
        dock, angle, cat = {}, {}, {}
        model_dict = {}
        for line in open(dlg):
            if line.startswith('DOCKED: MODEL'): # number of run
                model = int(line.split()[-1])

            elif 'Macromolecule file used to create Grid Maps' in line:
                pdbqt = str(re.split(r"\t", line)[-1]).strip()

            elif 'Coordinates of Central Grid Point of Maps' in line:
                center = eval(re.search(r"\(.*\)", line).group())

            elif line.startswith('DOCKED: ATOM'): # coordinates and energies of run
                number, atom, res, chain, x, y, z, vdw, Elec = line.replace('-', ' -').split()[2:11]

                ###  only args.residue atoms are considered
                if args.residue:
                    if res == args.residue:
                        dock.setdefault(model, []).append((float(vdw) + float(Elec)))
                else:
                    dock.setdefault(model, []).append((float(vdw) + float(Elec)))

                ###  only args.distance atom is considered   
                if args.distance:
                    if number == str(args.distance):
                        model_dict.setdefault(model, distance(center,tuple(map(float,(x,y,z)))))
                else:
                    model_dict.setdefault([model, 0])

                ###  only args.angle atom is considered   
                if args.angle:
                    if number in args.angle:
                        angle.setdefault(model, []).append(list(map(float,(x,y,z))))
                else:
                    angle.setdefault(model, 0)                

                ###  only args.catalytic atoms are considered    
                if args.catalytic:
                    if number in args.catalytic:
                        cat.setdefault(model, []).append(tuple(map(float,(x,y,z))))
                else:
                    cat.setdefault(model, [0])

            elif 'RANKING' in line: # parsing rmsd table
                tab.append(list(map(float,line.split()[:4])))

        df_dock = pd.DataFrame.from_dict(dock, orient='index')
        df_dock['vdw+Elec'] = df_dock.sum(axis=1) # sum of non binding energy for each run
        try:
            df_dist = pd.DataFrame.from_dict(model_dict, orient='index', columns = ['Distance']) # distance from coordinates (e.g. C1)
            df_cat = pd.DataFrame.from_dict(cat, orient='index', columns = args.catalytic) # coordinates of catalytic competence atoms
            df_angle = pd.DataFrame.from_dict(angle, orient='index', columns = args.angle) # coordinates for dihedral

            df_angle['dihedral'] = df_angle.apply(lambda v : abs(math.sin(calc_dihedral(Vector(v[0]), # coomputing sin(dihedral) from coordinates
                                                                      Vector(v[1]),
                                                                      Vector(v[2]),
                                                                      Vector(v[3])))), axis = 1)

            dfs = pd.concat([df_dock['vdw+Elec'], df_dist, df_cat, df_angle['dihedral']], axis=1) # concatenate dataframes
        except:
            dfs = pd.DataFrame(df_dock['vdw+Elec'])

        tab = pd.DataFrame(tab, columns = ['rank', 'subrank', 'run', 'energy']) # dataframe of HISTOGRAM in dlg 
        tab.set_index('run', inplace = True)

        df = pd.merge(tab, dfs, left_index=True, right_index = True)  # build total dataframe

        ### LARGEST CLUSTER ENERGY
        LC_sorted = df.groupby('rank').count().sort_values('subrank', ascending = False)['subrank'].index.astype(int).tolist()
        LC = df.groupby('rank').count()['subrank'].idxmax().astype(int)
        runs_in_LC = df.groupby('rank').count()['subrank'][LC_sorted[0]]
        LC_energies = (df[df['rank'] == LC_sorted[0]]['energy'],
                       df[df['rank'] == LC_sorted[0]]['vdw+Elec'])
        LCB = (LC_energies[0].min(), LC_energies[0].idxmin())
        LCM = (LC_energies[0].mean(), LC_energies[0].index.tolist())
        LCaaB = (LC_energies[1].min(), LC_energies[1].idxmin())
        LCaaM = (LC_energies[1].mean(), LC_energies[1].index.tolist())

        ### BEST CLUSTER ENERGY
        BC_sorted = df['rank'].drop_duplicates().astype(int).tolist()
        BC_energies = (df[df['rank'] == BC_sorted[0]]['energy'],
                       df[df['rank'] == BC_sorted[0]]['vdw+Elec'])
        BCB = (BC_energies[0].min(), BC_energies[0].idxmin())
        BCM = (BC_energies[0].mean(), BC_energies[0].index.tolist())
        BCaaB = (BC_energies[1].min(), BC_energies[1].idxmin())
        BCaaM = (BC_energies[1].mean(), BC_energies[1].index.tolist())

        try:
            runs_in_second_LC = df.groupby('rank').count()['subrank'][LC_sorted[1]]
        except:
            runs_in_second_LC = 0

        ratio_LC = round((1-(runs_in_second_LC/runs_in_LC))*100,2)

        try: ##try to open pdbqt
            pdb = PandasPdb().read_pdb(f"{dlg_dir}/{pdbqt}") # open pdb with biopandas
        except: ##if pdb is present
            print('Please re-install biopandas from\nconda install -c conda-forge biopandas')

        ### CATAYTIC COMPETENCE
        cat_comp = df[(df['dihedral'] >= 0.9)
              & (df['Distance'] <= 5)
              & (df['rank'] == LC) 
             ] # filter by catalyitic competences

        if args.catalytic:
            for k, c in cat.items():
                if k in cat_comp.index.tolist():
                    around_int = []
                    for center in c:
                        distances = pdb.distance(xyz=center, records=('ATOM'))
                        site = pdb.df['ATOM'][pdb.distance(xyz=center, records=('ATOM',)) <= 4].astype(str) # threshold of 3.5 A 
                        site['distances'] = distances
                        site['run'] = [k]*len(site)
                        site = site.sort_values('distances')        

                        intorno = site[(site['residue_name'].isin(['HIS',
                                                                  'ASP',
                                                                  'GLU',
                                                                  'TYR'
                                                                 ])) & 
                                      (site['atom_name'].isin(['OH', 
                                                               'ND1', 'NE2', 
                                                               'OE1', 'OE2',
                                                               'OD1', 'OD2',
                                                              ]))
                                     ][['run','record_name','residue_name','residue_number', 'atom_name', 'distances']]

                        around_int.append(intorno)
                    try:
                        around.append(pd.concat(around_int).drop_duplicates(['residue_name', 'residue_number']))
                        count.append(pd.concat(around_int).drop_duplicates('record_name'))
                    except:
                        pass
            try:
                cc = (pd.concat(around)
                      .groupby(['residue_name','residue_number', 'atom_name'])
                      .count()
                      .apply(list)
                      .sort_values('run', ascending = False))

                count = pd.concat(count).set_index('run')
                frequent_atom = ' '.join(list(zip(cc.index.tolist(), cc.record_name.tolist()))[0][0])
                frequent_atom_num = list(zip(cc.index.tolist(), cc.record_name.tolist()))[0][1]
            except:
                frequent_atom, frequent_atom_num = None, None
                count = pd.DataFrame([])

            ### CATAYTIC COMPETENCE
            cat_comp = cat_comp.merge(count, right_index = True, left_index = True, how = 'inner')
            cat_comp_num = len(cat_comp) ## number of catalytic runs in largest cluster

        else:
            frequent_atom, frequent_atom_num, cat_comp_num = 'None', 'None', 'None'

        print(os.path.basename(dlg), args.residue, 
              *BCB, *LCB, 
              *BCaaB, *LCaaB,
              *BCM, *LCM, 
              *BCaaM, *LCaaM,
              LC, runs_in_LC, ratio_LC, 
              df[df['rank'] == LC].mean().Distance, df[df['rank'] == LC].mean().dihedral,
              frequent_atom, frequent_atom_num, cat_comp_num,
              sep = '\t', 
              file=open(f"{dlg_dir}/{os.path.basename(dlg_dir)}.tsv", 'a'),
             )
    except:
        print(os.path.basename(dlg), sep = '\t', file=open(f"{dlg_dir}/{os.path.basename(dlg_dir)}.tsv", 'a'))

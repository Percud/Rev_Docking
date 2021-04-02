####    make table from dlg directory with different kind of energy    #### 
####    USAGE: ./get_atom_energy.py MY_RES_IN_LIGAND MY_DLG_DIRECTORY

import glob, re, sys, os, pandas as pd, argparse

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
                    help = 'Define FASTA input, gapped or not')
parser.add_argument('-r', '--residue',
                    help = 'Define residue three-letter code to compute energy separately')
parser.add_argument('-d', '--distance',
                    help = 'Define atom number to measure distance from grid center')
parser.add_argument('-c', '--catalytic',
                    nargs = '+',
                    help = 'Define atom number to measure distance from grid center')
args = parser.parse_args()

parser.add_argument('-o', '--output',
                    default = os.path.dirname(args.input[0]) + '/' + os.path.basename(os.path.dirname(args.input[0])) + '.tsv',
                    help = 'Define output filename of tsv format table')
args = parser.parse_args()


dlg_dir = os.path.dirname(args.input[0])  ## DLGs directory    


## print header of table
print('dlg file', 'Ligand', 
      'BCB (kcal/mol)', 'Run', 'LCB (kcal/mol)', 'Run', 
      'BCaaB (kcal/mol)', 'Run', 'LCaaB (kcal/mol)', 'Run', 
      'BCM (kcal/mol)', 'LCM (kcal/mol)', 
      'BCaaM (kcal/mol)', 'LCaaM (kcal/mol)',
      'LC', 'Num in LC', '1LC/2LC %', 'Distance (A)', 'Catalytic Energy',
      sep = '\t',
      file=open(args.output, 'w')) # tab headers

for dlg in args.input:
    try:
        tab, dock, model_dict, oh = [], {}, {}, {}
        for line in open(dlg):
            if line.startswith('DOCKED: MODEL'): # number of run
                model = int(line.split()[-1])

            elif 'Coordinates of Central Grid Point of Maps' in line:
                center = eval(re.search(r"\(.*\)", line).group())

            elif line.startswith('DOCKED: ATOM'): # coordinates and energies of run
                number, atom, res, chain, x, y, z, vdw, Elec = line.split()[2:11]
                
                ###  only args.residue atoms are considered
                if args.residue:
                    if res == args.residue:
                        dock.setdefault(model, []).append((float(vdw) + float(Elec)))
                else:
                    dock.setdefault(model, []).append((float(vdw) + float(Elec)))
                
                ###  only args.distance atom is considered   
                if args.distance:
                    if number == args.distance:
                        model_dict.setdefault(model, distance(center,tuple(map(float,(x,y,z)))))
                else:
                    model_dict.setdefault(model, 0)
                
                ###  only args.prova atoms are considered    
                if args.catalytic:
                    if number in args.catalytic:
                        oh.setdefault(model, []).append((float(vdw) + float(Elec)))
                else:
                    oh.setdefault(model, [0])

            elif 'RANKING' in line: # parsing rmsd table
                tab.append(list(map(float,line.split()[:4])))

        df_dock = pd.DataFrame.from_dict(dock, orient='index')
        df_dock['vdw+Elec'] = df_dock.sum(axis=1) # sum of non binding energy for each run

        df = pd.DataFrame(tab, columns = ['rank', 'subrank', 'run', 'energy'], dtype = int)
        df.set_index(['rank','run'], inplace = True)
        df = df.join(df_dock['vdw+Elec'], on='run').round(2)
        bo = df.groupby('rank')
        bo.idxmin().at[1,'vdw+Elec'][1]
        best_energy = (df.loc[1].iloc[0]['energy'], 
                       df.loc[1]['energy'].idxmin())
        largest_energy = (df.loc[df['subrank'].idxmax()[0]].iloc[0]['energy'], 
                          df.loc[df['subrank'].idxmax()[0]]['energy'].idxmin())
        best_res = (df.loc[1].iloc[0]['vdw+Elec'], 
                    df.loc[1]['subrank'].idxmin())
        largest_res = (df.loc[df['subrank'].idxmax()[0]].iloc[0]['vdw+Elec'], 
                       df.loc[df['subrank'].idxmax()[0]]['subrank'].idxmin())
        mean_largest = round(df.loc[df['subrank'].idxmax()[0]]['energy'].mean(),2)
        mean_best = round(df['energy'][1].mean(),2)
        best_largest_res = (df.loc[df['subrank'].idxmax()[0]]['vdw+Elec'].min(), 
                            df.loc[df['subrank'].idxmax()[0]]['vdw+Elec'].idxmin())
        mean_largest_res = round(df.loc[df['subrank'].idxmax()[0]]['vdw+Elec'].mean(),2)
        mean_best_res = round(df['vdw+Elec'][1].mean(),2)
        avg = round(df['vdw+Elec'].mean(),2)
        bestbest = (df['vdw+Elec'].min(), 
                    df['vdw+Elec'].idxmin()[1])
        largest_cluster = df['subrank'].idxmax()[0]
        largest_cluster_num = df['subrank'].max()
        try:
            second_largest_cluster_num = pd.DataFrame(bo.count().sort_values('subrank',ascending=False)).loc[2,'subrank']
        except:
            second_largest_cluster_num = 0
        ratio_num_LC=round((1-(second_largest_cluster_num/largest_cluster_num))*100,2)
        print(os.path.basename(dlg), args.residue, 
              *best_energy, *largest_energy, 
              *best_res, *largest_res, 
              mean_best, mean_largest,
              mean_best_res, mean_largest_res,
              largest_cluster, largest_cluster_num, ratio_num_LC, model_dict.get(largest_energy[1]), sum(oh.get(largest_energy[1])),
              sep = '\t', 
              file=open(args.output, 'a'),
             )
    except:
        print(os.path.basename(dlg), sep = '\t', file=open(args.output, 'a'))

    

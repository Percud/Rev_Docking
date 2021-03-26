####    make table from dlg directory with different kind of energy    #### 
####    USAGE: ./get_atom_energy.py MY_RES_IN_LIGAND MY_DLG_DIRECTORY

import glob, re, sys, os, pandas as pd

def distance(a,b):
    '''Return max distance given arrays in xyz of points'''
    x = [a[0], b[0]]
    y = [a[1], b[1]]
    z = [a[2], b[2]]
    return (((max(x)-min(x))**2+(max(y)-min(y))**2+(max(z)-min(z))**2)**0.5)

try:
    RES = str(sys.argv[1]) ## three-letter code of selected ligand atom
    dlg_dir = sys.argv[2]  ## DLGs directory
except:
    print('''
Returns a Dataframe.tsv given a three-letter code of the ligand and a dlg directory

USAGE : DLGdf three_letter_code dlg_directory

Example : DLGdf GLU docking_result
''')
    sys.exit(1)

    
cwd = os.getcwd()
os.chdir(dlg_dir)
dlg_dir = os.path.basename(os.getcwd())

## print header of table
print('dlg file', 'Ligand', 
      'BCB (kcal/mol)', 'Run', 'LCB (kcal/mol)', 'Run', 
      'BCaaB (kcal/mol)', 'Run', 'LCaaB (kcal/mol)', 'Run', 
      'BCM (kcal/mol)', 'LCM (kcal/mol)', 
      'BCaaM (kcal/mol)', 'LCaaM (kcal/mol)',
      'LC', 'Num in LC', '1LC/2LC %', 'Distance (A)'
      sep = '\t',
      file=open(dlg_dir+'.tsv', 'w')) # tab headers

for dlg in glob.glob('*.dlg'):
    try:
        tab, dock = [], {}
        model_dict = {}
        for line in open(dlg):
            if line.startswith('DOCKED: MODEL'): # number of run
                model = int(line.split()[-1])
            elif 'Coordinates of Central Grid Point of Maps' in line:
                center = eval(re.search(r"\(.*\)", line).group())
            elif line.startswith('DOCKED: ATOM'): # coordinates and energies of run
                number, atom, res, chain, x, y, z, vdw, Elec = line.split()[2:11]
                if res == RES:  # only RES atom are considered
                    dock.setdefault(model, []).append((float(vdw) + float(Elec)))
                if number == '2':
                    model_dict.setdefault(model, distance(center,tuple(map(float,(x,y,z)))))
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
        print(os.path.basename(dlg), RES, 
              *best_energy, *largest_energy, 
              *best_res, *largest_res, 
              mean_best, mean_largest,
              mean_best_res, mean_largest_res,
              largest_cluster, largest_cluster_num, ratio_num_LC, model_dict.get(largest_energy[1]),
              sep = '\t', 
              file=open(dlg_dir + '.tsv', 'a'),
             )
    except:
        print(os.path.basename(dlg), sep = '\t', file=open(dlg_dir + '.tsv', 'a'))

os.chdir(cwd)

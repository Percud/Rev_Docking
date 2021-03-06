#!/usr/bin/env python3

####    make table from dlg directory with different kind of energy    #### 
####    USAGE: ./get_atom_energy.py MY_RES_IN_LIGAND MY_DLG_DIRECTORY

import glob, re, sys, os, pandas as pd

try:
    RES = str(sys.argv[1]) ## three-letter code of selected ligand atom
    dlg_dir = sys.argv[2]
except:
    print('Returns a Dataframe.tsv given a three-letter code of the ligand and a dlg directory',
          '\nUSAGE : DLGdf three_letter_code dlg_directory',
          '\nExample : DLGdf GLU docking_result')
    sys.exit(1)
    
## print header of table
print('dlg file', 'Ligand', 
      'BCB (kcal/mol)', 'Run', 'LCB (kcal/mol)', 'Run', 
      'BCaaB (kcal/mol)', 'Run', 'LCaaB (kcal/mol)', 'Run', 
      'BCM (kcal/mol)', 'LCM (kcal/mol)', 
      'BCaaM (kcal/mol)', 'LCaaM (kcal/mol)',
      'LC', 'Num in LC', '1LC/2LC %',
      sep = '\t',
      file=open(dlg_dir+dlg_dir.split('/')[-2]+'.tsv', 'w')) # tab headers

for dlg in glob.glob(dlg_dir+'/*.dlg'):
    try:
        tab, dock = [], {}
        for line in open(dlg):
            if line.startswith('DOCKED: MODEL'):
                model = int(line.split()[-1])
            elif line.startswith('DOCKED: ATOM'):
                number, atom, res, chain, x, y, z, vdw, Elec = line.split()[2:11]
                if res == RES:
                    dock.setdefault(model, []).append((float(vdw)+float(Elec)))
            elif 'RANKING' in line:
                tab.append(list(map(float,line.split()[:4])))

        df_dock = pd.DataFrame.from_dict(dock, orient='index')
        df_dock['vdw+Elec'] = df_dock.sum(axis=1)

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
              largest_cluster, largest_cluster_num, ratio_num_LC,
              sep = '\t', file=open(dlg_dir+dlg_dir.split('/')[-2]+'.tsv', 'a'))
    except:
        print(os.path.basename(dlg), sep = '\t', file=open(dlg_dir+dlg_dir.split('/')[-2]+'.tsv', 'a'))

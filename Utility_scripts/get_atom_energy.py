#!/usr/bin/env python3

####    make table from dlg directory with different kind of energy    #### 
####    USAGE: ./get_atom_energy.py MY_RES_IN_LIGAND MY_DLG_DIRECTORY

import glob, re, sys, os, pandas as pd

RES = str(sys.argv[1]) ## three-letter code of selected ligand atom
dlg_dir = sys.argv[2]

## print header of table
print('dlg file', 'Ligand', 'BC (kcal/mol)', 'Run', 'LC (kcal/mol)', 'Run', 'BCaa (kcal/mol)', 'Run', 'LCaa (kcal/mol)', 'Run', 'LCaaB (kcal/mol)', 'Run', 'LCaaM (kcal/mol)', 'BB (kcal/mol)', 'Run', 'LC', 'Num in LC', sep = '\t', file=open(dlg_dir+dlg_dir.split('/')[-2]+'.tsv', 'w')) # tab headers
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
        best_energy = (df.loc[1].iloc[0]['energy'], df.loc[1]['energy'].idxmin())
        largest_energy = (df.loc[df['subrank'].idxmax()[0]].iloc[0]['energy'], df.loc[df['subrank'].idxmax()[0]]['energy'].idxmin())
        best_res = (df.loc[1].iloc[0]['vdw+Elec'], df.loc[1]['subrank'].idxmin())
        largest_res = (df.loc[df['subrank'].idxmax()[0]].iloc[0]['vdw+Elec'], df.loc[df['subrank'].idxmax()[0]]['subrank'].idxmin())
        best_largest_res = (df.loc[df['subrank'].idxmax()[0]]['vdw+Elec'].min(), df.loc[df['subrank'].idxmax()[0]]['vdw+Elec'].idxmin())
        mean_largest_res = round(df.loc[df['subrank'].idxmax()[0]]['vdw+Elec'].mean(),2)
        avg = round(df['vdw+Elec'].mean(),2)
        bestbest = (df['vdw+Elec'].min(), df['vdw+Elec'].idxmin()[1])
        largest_cluster = df['subrank'].idxmax()[0]
        largest_cluster_num = df['subrank'].max()
        print(os.path.basename(dlg), RES, *best_energy, *largest_energy, *best_res, *largest_res, *best_largest_res, mean_largest_res, *bestbest, largest_cluster, largest_cluster_num, sep = '\t', file=open(dlg_dir+dlg_dir.split('/')[-2]+'.tsv', 'a'))
    except:
        print(os.path.basename(dlg), sep = '\t', file=open(dlg_dir+dlg_dir.split('/')[-2]+'.tsv', 'a'))
        



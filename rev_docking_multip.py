#!/usr/bin/env python3

import re, sys, os, subprocess, glob, configparser
from multiprocessing import Pool

cwd = os.getcwd()

#####    config parser    ####
config = configparser.ConfigParser(inline_comment_prefixes="#")
config.optionxform = str
config.read(sys.argv[1])
locals().update(dict(config.items('DEFAULT')))

####    call to scripts    ####
pythonsh = MGL_ROOT+'/bin/pythonsh '
u24 = MGL_ROOT+'/MGLToolsPckgs/AutoDockTools/Utilities24/'
def prepare_receptor(inf,outf):
    os.system(pythonsh+u24+'prepare_receptor4.py -A "checkhydrogens" -e -r '+inf+' -o '+outf)
def prepare_gpf(receptor, ligand, output):
    os.chdir(out_dir)
    os.system(pythonsh+u24+'prepare_gpf4.py -i reference.gpf -r '+receptor+' -l '+ligand+' -o '+output)
    os.chdir(cwd)
def prepare_dpf(receptor, ligand, output):    
    os.chdir(out_dir)
    os.system(pythonsh+u24+'prepare_dpf42.py -i reference.dpf -r '+receptor+' -l '+ligand+' -o '+output)
    os.chdir(cwd)
def autogrid(gpf, glg):
    os.chdir(out_dir)
    os.system('autogrid4 -p '+gpf+' -l '+glg)
    os.chdir(cwd)
def autodock(dpf, dlg):
    os.chdir(out_dir)
    os.system('autodock4 -p '+dpf+' -l '+dlg)
    os.chdir(cwd)

####   make outdir   ####
i = 1
while os.path.exists(outdir+'_%s' % i):
    i += 1
out_dir = outdir+'_%s' % i+'/'
os.makedirs(out_dir)
if os.path.exists(outdir+'_%s' % i):
    print('writing output in '+out_dir)

####    read coordinates    ####
coordinates = {}
try:
    for line in open(coord_file):
        if '#' in line: #skip comments
            next
        else:
            name, x, y, z = line.split()
            coordinates.setdefault(name, []).append([x,y,z])
    print('Reading '+coord_file+' ...')
except:
    print(coord_file+" doesn't exist")

try:
    os.system('cp '+ligand+' '+out_dir)
    lig_name = os.path.splitext(os.path.basename(ligand))[0]
except:
    print(ligand+" doesn't exist")


#for pdb, xyz in coordinates.items():
#    pdb_name = os.path.splitext(pdb)[0]
#    for c in xyz:
#        ####    prepare receptor, gpf and dpf  ####
#        os.system('cp '+gpf_par+' '+dpf_par+' '+out_dir)
#
#        gpf = out_dir+os.path.basename(gpf_par) # define new path
#        dpf = out_dir+os.path.basename(dpf_par)
#
#        with open(gpf, 'a') as gpf, open(dpf, 'a') as dpf: 
#            gpf.write('\ngridcenter '+' '.join(map(str,c))) # write gridcenter
#        i=1
#        while os.path.exists(out_dir+pdb_name+'_%s.pdbqt' %i):
#            i += 1 # incrementing file number
#
#        prepare_receptor(receptor_dir+'/'+pdb, out_dir+pdb_name+'_%s.pdbqt' %i)
#        prepare_gpf(pdb_name+'_%s.pdbqt' %i, os.path.basename(ligand), pdb_name+'_%s.gpf' % i)
#        prepare_dpf(pdb_name+'_%s.pdbqt' %i, os.path.basename(ligand), pdb_name+'_'+lig_name+'_%s.dpf' % i)
#
#        ## solve bug space in dpf ##
#        print(re.sub('#', ' #', open(out_dir+pdb_name+'_'+lig_name+'_%s.dpf' % i).read()), file=open(out_dir+pdb_name+'_'+lig_name+'_%s.dpf' % i, 'w'))
#
#        ####    autogrid    ####
#        autogrid(pdb_name+'_%s.gpf' % i, pdb_name+'_%s.glg' % i)
#        
#        pid = os.fork()
#        if pid == 0:
#            ####    autodock    ####
#            autodock(pdb_name+'_'+lig_name+'_%s.dpf' % i, pdb_name+'_'+lig_name+'_%s.dlg' % i)
#

def reverse_docking(pdb):  
    xyz = coordinates.get(pdb)
    pdb_name = os.path.splitext(pdb)[0]
    for c in xyz:
        
        ####    prepare receptor, gpf and dpf  ####
        os.system('cp '+gpf_par+' '+dpf_par+' '+out_dir)

        gpf = out_dir+os.path.basename(gpf_par) # define new path
        dpf = out_dir+os.path.basename(dpf_par)

        with open(gpf, 'a') as gpf, open(dpf, 'a') as dpf: 
            gpf.write('\ngridcenter '+' '.join(map(str,c))) # write gridcenter
            dpf.write('\ngridcenter '+' '.join(map(str,c)))
        i=1
        while os.path.exists(out_dir+pdb_name+'_%s.pdbqt' %i):
            i += 1 # incrementing file number

        print(open(receptor_dir+'/'+pdb, 'r').read().replace('HETATM', '#HETATM'), file=open(out_dir+pdb, 'w')) # comment HETATM
        prepare_receptor(receptor_dir+'/'+pdb, out_dir+pdb_name+'_%s.pdbqt' %i)
        prepare_gpf(pdb_name+'_%s.pdbqt' %i, os.path.basename(ligand), pdb_name+'_%s.gpf' % i)
        prepare_dpf(pdb_name+'_%s.pdbqt' %i, os.path.basename(ligand), pdb_name+'_'+lig_name+'_%s.dpf' % i)

        ## solve bug space in dpf ##
        print(re.sub('#', ' #', open(out_dir+pdb_name+'_'+lig_name+'_%s.dpf' % i).read()), file=open(out_dir+pdb_name+'_'+lig_name+'_%s.dpf' % i, 'w'))

        ####    autogrid    ####
#        autogrid(pdb_name+'_%s.gpf' % i, pdb_name+'_%s.glg' % i)

        ####    autodock    ####
#        autodock(pdb_name+'_'+lig_name+'_%s.dpf' % i, pdb_name+'_'+lig_name+'_%s.dlg' % i)


if __name__ == '__main__':
    pool = Pool(20)                       # Create a multiprocessing Pool
    pool.map(reverse_docking,coordinates)  # process data_inputs iterable with pool
    



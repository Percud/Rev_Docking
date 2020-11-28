#!/home/marco/anaconda3/bin/python3
from revdock import *

cwd=os.getcwd()

# ## HUMAN PLPome          
# ncbi_acc=pd.read_csv('http://bioinformatics.unipr.it/B6db/tmp/Homo_sapiens.tab',sep='\t',header=None)[2].tolist()
# accession=convert_ac(ncbi_acc,'P_REFSEQ_AC','ACC').To.tolist()
# get_models(accession,'9606','Human_PLP_swissmodel')

# ## Download fasta from uniprot
# os.chdir('Human_PLP_swissmodel_9606/')
# for a in set(accession):
#     try:
#         wget.download('https://www.uniprot.org/uniprot/'+a+'.fasta')
#     except:
#         pass
    
# ## Get coord from catalytic lysine
# output=[]
# for pdb in glob.glob('*.pdb'):
#     fa = pdb.strip('.pdb')+'.fa'
#     id = os.path.basename(fa).split('_')[0]
#     pdb2fasta(pdb, fa)
#     structure = PDB().get_structure(pdb, pdb)[0]
#     try:
#         features = pd.DataFrame(getfeatures(id).features[0])
#         lys = features[(features.category == 'PTM')&(features.description.str.contains('pyridoxal'))].begin.dropna().tolist()[0]
#         for k,v in convert(id+'.fasta', fa, lys).items():
#             coord=tuple(structure[k][v]['NZ'].get_coord())
#             output.append([pdb,k,v,*coord])
#     except:
#         output.append([pdb])
#     os.remove(fa)
#     os.remove(id+'.fasta')
# pd.DataFrame(output, columns=['pdb','chain','res','x','y','z']).to_csv('Human_coord.csv',index=False)
# os.chdir(cwd)

## MOUSE PLPome          
ncbi_acc=pd.read_csv('http://bioinformatics.unipr.it/B6db/tmp/Mus_musculus.tab',sep='\t',header=None)[2].tolist()
accession=convert_ac(ncbi_acc,'P_REFSEQ_AC','ACC').To.tolist()
get_models(accession,'10090','Mouse_PLP_swissmodel')

## Get coord from catalytic lysine
os.chdir('Mouse_PLP_swissmodel_10090/')
output=[]
for pdb in glob.glob('*.pdb'):
    ## LLP residue ##
    print(re.sub(r'HETATM(.*\s*[N|C|O|CA|CB|CG|CE|CD|NZ]\s*)LLP',
                 r'ATOM  \1LYS',open(pdb).read()),
          file=open(pdb,'w'))
    print(re.sub(r'(?m).*LLP.*\n?', '', open(pdb).read()),
          file=open(pdb,'w'))
    
    fa = pdb.split('.pdb')[0]+'.fa'
    id = os.path.basename(fa).split('_')[0]
    pdb2fasta(pdb, fa)
    structure = PDB().get_structure(pdb, pdb)[0]
    try:
        uni =  getfeatures(id)
        features = pd.DataFrame(uni.features[0])
        print('>'+id+'\n'+uni['sequence.sequence'].values[0], file=open(id+'.fasta', 'w'))
        lys = features[(features.category == 'PTM')&(features.description.str.contains('pyridoxal'))].begin.dropna().tolist()[0]
        for k,v in convert(id+'.fasta', fa, lys).items():
            coord=tuple(structure[k][v]['NZ'].get_coord())
            output.append([pdb,k,v,*coord])
    except:
        output.append([pdb])
    os.remove(fa)
    os.remove(id+'.fasta')

pd.DataFrame(output, columns=['pdb','chain','res','x','y','z']).to_csv('Mouse_coord.csv',index=False)
os.chdir(cwd)
import pandas as pd, tarfile, urllib, wget, os, glob, requests, numpy as np, re
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast import NCBIXML 
from Bio import SeqIO, SearchIO
from Bio.PDB import PDBParser as PDB, PDBIO

def convert(query, subject, num):
    """return number of subject fasta given a query number in dict format"""
    num=int(num)
    blastp(query=query, subject=subject, out=query+'_'+subject+'.xml', outfmt=5, max_hsps=1)()
    blastp(query=query, subject=subject, out=query+'_'+subject+'.txt', max_hsps=1)()
    xml = SearchIO.read(query+'_'+subject+'.xml', "blast-xml")
    d={}
    for n in range(len(xml)):
        x=xml[n][0]
        hit=list('X'*x.hit_start)+list(x.hit)
        q=list('X'*x.query_start)+list(x.query)
        hit_gap = np.array([index for index, value in enumerate(hit) if value == '-'])
        q_gap = np.array([index for index, value in enumerate(q) if value == '-'])
        hit[:] = (value for value in hit if value != '-')
        q[:] = (value for value in q if value != '-')
        sub=(num -len(hit_gap[hit_gap<num-x.query_start+x.hit_start])+len(q_gap[q_gap<num])-x.query_start+x.hit_start)
        name=(xml[n].description)
        if num-1 in hit_gap+x.query_start-x.hit_start:
            sub=None
            d.setdefault(name,{num:q[num-1],sub:sub})
        else:
            d.setdefault(name,{num:q[num-1],sub:hit[sub-1]})
    os.remove(query+'_'+subject+'.xml')
    return(d)

def convert_ac(list,From,to):
    """convert accession list"""
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': From,
    'to': to,
    'format': 'tab',
    'query': ' '.join(list)
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    return pd.read_csv(urllib.request.urlopen(req), sep='\t')

def pdb2fasta(pdb,fasta):
    """convert pdb file to fasta format for each chain"""
    with open(fasta,'w') as outfa:
        for record in SeqIO.parse(pdb, 'pdb-atom'):
            if record.annotations['start'] < 1:
                print('>' +record.annotations['chain']+'\n'+record.seq[1:],file = outfa)
            else:
                print('>' +record.annotations['chain']+'\n'+'X'*(record.annotations['start']-1)+record.seq,file = outfa)

def getfeatures(uniprot_id):
    """return the Uniprot complete dataset for a given Uniprot ID"""
    try:
        r = requests.get("https://www.ebi.ac.uk/proteins/api/proteins?size=-1&accession=" + uniprot_id, 
                         headers={ "Accept" : "application/json"})
        data = pd.json_normalize(r.json())
        return data
    except:
        return str(uniprot_id) + ' not_found'

def get_models(list,taxid,out):
    cwd = os.getcwd()
    os.mkdir(out+'_'+taxid)
    os.chdir(out+'_'+taxid)
    tar_swissmodel = urllib.request.urlopen('https://swissmodel.expasy.org/repository/download/core_species/'+taxid+'_meta.tar.gz')
    tar = tarfile.open(fileobj=tar_swissmodel, mode="r|gz")
    index = tar.extractall()
    df=pd.read_csv('SWISS-MODEL_Repository/INDEX',sep='\t',comment='#')
    tar.close()
    swissmodel=df[df.UniProtKB_ac.str.contains('|'.join(list))].drop_duplicates('url')
    swissmodel.to_csv(out+'.csv')
    for url in swissmodel.url.tolist():
        try:
            wget.download(url)
        except:
            with open('downlod.log','a') as log:
                print(url, 'Not found',file=log)
    os.chdir(cwd)
    
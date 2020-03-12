#!/bin/bash
#Replace HETATM for Lys atoms of the internal aldimine 

#Usage procedure_transform_LLP_aldimmine.sh 

for f in *.pdb
do
    #Replace LLP with LYS for lysine residue
    sed -i  's/HETATM\(.*\s*[N|C|O|CA|CB|CG|CE|CD|NZ]\s*\)LLP/ATOM  \1LYS/' $f
done


echo "Replaced LLP HETATM with ATOM in Lys internal aldimmine"

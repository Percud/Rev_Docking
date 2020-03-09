#!/bin/bash
#Replace HETATM for Lys atoms of the internal aldimine 

for f in $@
do
  sed -i  's/HETATM\(.*\s*[N|C|O|CA|CB|CG|CE|CD|NZ]\s*\)LLP/ATOM  \1LYS/' $f
done

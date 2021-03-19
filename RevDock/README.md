# Reverse Docking virtual screening of enzymes active sites.

Python scripts and data to performs a virtual screening of a library of enzyme active sites using a given substrate.
Coinceived for a High Performance Computer (HPC) cluster

## Getting Started


### Prerequisites

AutoDock, AutoGrid >=4.2.6

Mgltools >=1.5.6

Python3 >=3.6

## Example
Screen a library of enzyme active sites with a particular substrate
```
python3 reverse_docking_p.py rd_conf_HTML_PLP.txt
```

See  [rd_conf_HTML_PLP.txt](https://github.com/Percud/Rev_Docking/edit/master/rd_conf_HTML_PLP.txt) for required files

______________________________________________________________
A list of pdb entries with x y z coordinates to center the grid 
can be obtained with the following utility:

```
./Utility_scripts/procedure_get_coord_from_resn.pl -a NZ Human/Human_Orf_pred.resn Human/pdb
```
  The command will extract coordinates from the pdb files stored in the Huma/pdb directory.

  Pdb_id, chain, and residue number are specified in the input file; the -a option specifies the particular atom (default CA).

  See  [Human_Orf_pred.resn](https://github.com/Percud/Rev_Docking/edit/master/Human/Human_Orf_pred.resn) for input file specs

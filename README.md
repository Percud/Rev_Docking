# Project Title

Enzyme substrate active site reverse docking

## Getting Started


### Prerequisites

AutoDock >=4.2.6
Perl, Phyton, linux

## Example

```
./procedure_reverse_docking_p.pl rd_conf_HTML_PLP.txt
```

Lists of pdb entries with x y z coordinates to center the grid 
can be obtained with the following utility:

```
Utility_scripts/procedure_get_coord_from_resn.pl -a NZ Human/Human_Orf_pred.resn Human/pdb
```
The command will extract coordinates from the pdb files stored in the Huma/pdb directory
See  [Human_Orf_pred.resn](https://github.com/Percud/Rev_Docking/edit/master/Human/Human_Orf_pred.resn) for input file specs

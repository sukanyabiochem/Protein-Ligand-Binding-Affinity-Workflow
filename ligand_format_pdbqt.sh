#!/bin/bash

##This generates:
## ligand_1.pdb
##ligand_2.pdb
##ligand_3.pdb
##...

for num in {1..34}
do
mkdir -p ligand_prep ; ## create folder to make library of the ligand ## 
obabel ligand_${num}.pdb -O ligand_prep/ligand_${num}.mol2 -h ; ## convert the pdb to mol2 file ##  
obabel ../ligand_${num}.pdb -O ../ligand_${num}.sdf ## convert the pdb to sdf file format ##
mk_prepare_ligand.py -i ligand_prep/ligand_${num}.mol2 -o ligand_prep/ligand_${num}.pdbqt ; ## Convertion of the ligand format to autodoc compatable pdbqt ##  
echo "$num"
done

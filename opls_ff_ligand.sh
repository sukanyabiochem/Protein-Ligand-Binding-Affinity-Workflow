## automatic python command line to download ITP file of OPLS force-field from Ligpargen Server ##

for i in {1..34}
do	
python3 final_python/ligpargen.py -i ligand_${i}.pdb -o ligpergen_${i} --only prm,rtf,itp,gro,zip
done


# Download  OPLS force field for MD simulation run in GROMACS for one ligand ##
python ligpargen.py -i ligand_1.pdb -o ligpargen_output

# Only specific file types
python ligpargen.py -i ligand_1.pdb -o ligpargen_output --only prm,rtf,itp,gro,zip

# Charged molecule (cm1a + net charge)
python ligpargen.py -i ligand_1.pdb --charge-model cm1a --charge -1

# List download targets without saving
python ligpargen.py -i ligand_1.pdb --no-download

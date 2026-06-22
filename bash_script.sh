for num in {1..34}
do
obabel ligand_${num}.pdb -O ligand_${num}.mol2 -h       	
mk_prepare_ligand.py -i ligand_${num}.mol2 -o ligand_${num}.pdbqt
echo "$num"
done

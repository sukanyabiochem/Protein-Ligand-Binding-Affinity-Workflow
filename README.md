# Protein-Ligand-Interaction-Energy-ML


## create a conda enviorment for Autodoc ##
conda create -n vina python=3.10
conda activate vina
conda config --env --add channels conda-forge
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
conda install -c conda-forge vina

#check installation
vina --help
#two other package for file conversion from pdb to pdbqt
conda install -c conda-forge openbabel
conda install -c conda-forge meeko


##install ImageMagick for convert command and Convert PNG → GIF
#sudo apt update
sudo apt install imagemagick
convert -delay 5 -loop 0 frames/frame*.png animation.gif
convert -delay 10 frames/frame*.png animation.gif

##AUTODOCK
#split the 34 ligand from one file into different##
obabel 34_ligand_convert.pdb -O ligand_.pdb -m

## check their str in VMD ##
vmd ligand_1.pdb

## prepare the receptor file (Protein) ##
##mk_prepare_ligand.py -i ligand_1.pdb -o ligand_1.pdbqt

## convert the pdb to mol2 file ##
obabel ligand_1.pdb -O ligand_1.mol2 -h

vmd ligand_1.mol2

## Convert the PDB --> mol2 & prepare the pdbqt file for the all 34 ligands ## 
##vi bash_script.sh 
bash bash_script.sh

## create the directory ##
mkdir output

## Prepare the pdbqt file for the Protein ##
mk_prepare_receptor.py -i FVIIa.pdb -o FVIIa_receptor -p -v -g --box_size 50 50 50 --box_center 23.768259048461914 16.592430114746094 12.018509864807129

## creating the map file ##
autogrid4 -p FVIIa_receptor.gpf -l FVIIa_receptor.glg

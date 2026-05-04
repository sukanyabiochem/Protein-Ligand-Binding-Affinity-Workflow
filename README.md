# Protein-Ligand-Interaction-Energy-ML


## create a conda enviorment for Autodoc ##
1. conda create -n vina python=3.10
2. conda activate vina
3. conda config --env --add channels conda-forge
4. conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
5. conda install -c conda-forge vina

#check installation
vina --help
#two other package for file conversion from pdb to pdbqt
1. conda install -c conda-forge openbabel
2. conda install -c conda-forge meeko


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

## Running Atuodoc using Autodock4 forcefield ##

vina  --ligand ligand_1.pdbqt --maps FVIIa_receptor --scoring ad4 --exhaustiveness 32 --out output/FVIIa_ligand_ad4_out.pdbqt

##movie in pymol

##create folder in pymol using command:
1. import os
2. os.makedirs("frames", exist_ok=True)

#you already have a trajectory loaded, you can just use frames:
1. mplay
2. png frame, width=800, height=600, dpi=300, ray=1
3. mpng frames/frame

##better quality
1. set ray_trace_frames, 1
2. mpng frames/frame


##install ImageMagick for convert command and Convert PNG → GIF
#sudo apt update
1. sudo apt install imagemagick
2. convert -delay 5 -loop 0 frames/frame*.png animation.gif
3. convert -delay 10 frames/frame*.png animation.gif

## combine the FVIIa and the drug molecules files ##

# convert the download GRO files of ligand to PDB format #

##for i in {4..10}
##do
##grofile=$(find ligand/ligpargen_${i} -maxdepth 1 -type f -name "*.gro" 2>/dev/null | head -n 1)
##outfile="ligand/ligpargen_${i}/ligand_${i}.pdb"
##	echo "Converting $grofile -> $outfile"
##        gmx editconf -f "$grofile" -o "$outfile"
##done

for i in 4
do	
mkdir -p P_ligand_${i}
cp -r P_ligand_2/input-file P_ligand_${i}
cp sub1.sbatch P_ligand_${i} 
cp Protein/topol.top Protein/posre.itp P_ligand_${i}
gmx editconf -f Protein/FVIIa_updated.gro -o Protein/FVIIa_updated.pdb
cat Protein/FVIIa_updated.pdb ligand/ligpargen_${i}/ligand_${i}.pdb > P_ligand_${i}/Pro_Lig_${i}.pdb

cat > P_ligand_${i}/ions.mdp <<EOF
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF
## automatic job submssion of including minimization, equilibration, production run etc ##
cat > P_ligand_${i}/sub1.sh <<EOF

    # Step 1: Energy minimization
gmx grompp -f input-file/minim.mdp -c P_Lig_${i}_ions.gro -r P_Lig_${i}_ions.gro -p topol.top -n group.ndx -o equil1.tpr -maxwarn 1 
gmx mdrun -deffnm equil1 -ntmpi 30

gmx grompp -f input-file/nvt.mdp -c equil1.gro -r equil1.gro -p topol.top -n group.ndx -o equil2.tpr -maxwarn 1
gmx mdrun -deffnm equil2 -ntmpi 1

gmx grompp -f input-file/npt.mdp -c equil2.gro -r equil2.gro -t equil2.cpt -p topol.top -n group.ndx -o equil3.tpr -maxwarn 1
gmx mdrun -deffnm equil3 -ntmpi 1

gmx grompp -f input-file/equil4.mdp -c equil3.gro  -t equil3.cpt -p topol.top -n group.ndx -o equil4.tpr
gmx mdrun -deffnm equil4 -ntmpi 1

gmx grompp -f input-file/equil5.mdp -c equil4.gro -t equil4.cpt -p topol.top -n group.ndx -o equil5.tpr
gmx mdrun -deffnm equil5 -ntmpi 1

gmx grompp -f input-file/equil6.mdp -c equil5.gro -t equil5.cpt -p topol.top -n group.ndx -o equil6.tpr
gmx mdrun -deffnm equil6 -ntmpi 1
EOF


cd P_ligand_${i}
sed -i '/END/d' Pro_Lig_${i}.pdb
sed -i '/CRYST1/d' Pro_Lig_${i}.pdb
gmx editconf -f Pro_Lig_${i}.pdb -resnr 1 -o Pro_Lig_${i}.gro
gmx editconf -f Pro_Lig_${i}.gro -o boxed_P_Lig_${i}.gro -c -d 1.2 -bt cubic



ligfolder="ligand/ligpargen_${i}"
itpfile=$(basename "$ligfolder"/ligand_${i}.itp)
##itpfile=$(find "$ligfolder" -maxdepth 1 -type f -name "*.itp" -printf "%f\n" | head -n 1)
itpfilef=$(basename "$ligfolder"/forcefield_ligand${i}.itp)

sed -i "/#endif/a\\
; Include forcefield parameters ligand\\
#include \"../$ligfolder/$itpfile\"
" topol.top

sed -i '/Protein[[:space:]]*1/a\
UNL                 1
' topol.top

##sed -i "/#include "oplsaa.ff\/forcefield.itp"/a\
##include "../$ligfolder/$itpfilef\"
##" topol.top
##sed -i '/#include "oplsaa.ff\/forcefield.itp"/a\
##include "../'"$ligfolder"'/'"$itpfilef"'"
##' topol.top
sed -i '/#include "oplsaa.ff\/forcefield.itp"/a\
#include "../'"$ligfolder"'/'"$itpfilef"'"
' topol.top

## add solvant ##
gmx solvate -cp boxed_P_Lig_${i}.gro -cs spc216.gro -o P_Lig_${i}_solvated.gro -p topol.top

sed -i '/ligand_.*\.itp/a\
#include "oplsaa.ff/tip3p.itp"
' topol.top

## add ions ##
gmx grompp -f ions.mdp -c P_Lig_${i}_solvated.gro -p topol.top -o P_Lig_${i}_ions.tpr

echo -e "15" | gmx genion -s P_Lig_${i}_ions.tpr -o P_Lig_${i}_ions.gro -p topol.top -pname NA -nname CL -neutral

sed -i '/#include "oplsaa.ff\/tip3p.itp"/a\
#include "oplsaa.ff/ions.itp"
' topol.top



## generate index file ##
echo -e "1|13 \n q" | gmx make_ndx -f P_Lig_${i}_ions.gro -o group.ndx 


## automatic job submssion of including minimization, equilibration, production run etc ##
bash sub1.sh
cd ../


done

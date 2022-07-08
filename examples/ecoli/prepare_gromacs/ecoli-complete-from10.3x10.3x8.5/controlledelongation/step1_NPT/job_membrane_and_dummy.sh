#!/bin/bash
#SBATCH --job-name=membrane_and_dummy
#SBATCH --account=nn4654k
#SBATCH --time=0-12:00:0
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128

# Safety options, they stop the script if any error occours.
set -o errexit  # Recommended for easier debugging
set -o nounset  # Treat unset variables as errors

# Load modules
module purge
module load GROMACS/2021.5-foss-2021b

# export paths of files that do not need changing
INPUT_PATH="/cluster/projects/nn4654k/samiran/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-complete-from10.3x10.3x8.5"
# modify topol.top
#cp ${INPUT_PATH}/topol.top .
#echo '#include "toppar/martini_v2.0_lipids_all_201506.itp"' | cat - topol.top > temp && mv temp topol.top
#echo '#include "toppar/martini_dummy.itp"'              | cat - topol.top > temp && mv temp topol.top

cp -r ${INPUT_PATH}/{index.ndx,ecoli_ionised.gro} topol.top toppar mdp ${SCRATCH}/.
cd $SCRATCH
mkdir -p out

#run
# Minimization
set env GMX_MAXCONSTRWARN -1
# step6.0: softcore minimization
gmx grompp -f mdp/step6.0_minimization.mdp -p topol.top -n index.ndx -c ecoli_ionised.gro -o out/step6.0_minimization.tpr
gmx mdrun -deffnm out/step6.0_minimization #-v

# step6.1: no softcore minimization
gmx grompp -f mdp/step6.1_minimization.mdp -p topol.top -n index.ndx -c out/step6.0_minimization.gro -o out/step6.1_minimization.tpr
gmx mdrun -deffnm out/step6.1_minimization #-v

unset env GMX_MAXCONSTRWARN

# Equilibration
# step 6.2: equilibration NVT timestep = 0.002 ps
gmx grompp -f mdp/step6.2_equilibration.mdp -p topol.top -n index.ndx -c out/step6.1_minimization.gro -r out/step6.1_minimization.gro -o out/step6.2_equilibration.tpr
gmx mdrun -deffnm out/step6.2_equilibration -rdd 1.3 -v


## step 6.3: equilibration NVT timestep = 0.005 ps
#gmx grompp -f mdp/step6.3_equilibration.mdp -p topol.top -n index.ndx -c out/step6.2_equilibration.gro -r out/step6.2_equilibration.gro -o out/step6.3_equilibration.tpr
#gmx mdrun -deffnm out/step6.3_equilibration -rdd 1.3 -v
#
## step 6.4: equilibration NVT timestep = 0.010 ps
#gmx grompp -f mdp/step6.4_equilibration.mdp -p topol.top -n index.ndx -c out/step6.3_equilibration.gro -r out/step6.3_equilibration.gro -o out/step6.4_equilibration.tpr
#gmx mdrun -deffnm out/step6.4_equilibration -rdd 1.3 -v
#
## step 6.5: equilibration NVT timestep = 0.015 ps
#gmx grompp -f mdp/step6.5_equilibration.mdp -p topol.top -n index.ndx -c out/step6.4_equilibration.gro -r out/step6.4_equilibration.gro -o out/step6.5_equilibration.tpr
#gmx mdrun -deffnm out/step6.5_equilibration -rdd 1.3 -v
#
## step 6.6: equilibration NVT timestep = 0.020 ps
#gmx grompp -f mdp/step6.6_equilibration.mdp -p topol.top -n index.ndx -c out/step6.5_equilibration.gro -r out/step6.5_equilibration.gro -o out/step6.6_equilibration.tpr
#gmx mdrun -deffnm out/step6.6_equilibration -rdd 1.3 -v


# Copy back results from $SCRATCH to submission directory.
cp -r out $SLURM_SUBMIT_DIR/gromacs_out

exit 0

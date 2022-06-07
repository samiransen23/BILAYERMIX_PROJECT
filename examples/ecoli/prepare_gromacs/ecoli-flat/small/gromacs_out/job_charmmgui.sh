#!/bin/bash
#SBATCH --job-name=pope-popg-ca
#SBATCH --account=nn4654k
#SBATCH --time=0-0:30:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --qos=devel

# Safety options, they stop the script if any error occours.
set -o errexit  # Recommended for easier debugging
set -o nounset  # Treat unset variables as errors

# Load modules
module purge
module load GROMACS/2021.5-foss-2021b

cp -r gromacs $SCRATCH/.
cd $SCRATCH/gromacs

#run
# Minimization
set env GMX_MAXCONSTRWARN -1
# step6.0 - soft-core minimization
gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c step5_charmm2gmx.pdb -r step5_charmm2gmx.pdb -p system.top -n index.ndx
gmx mdrun -deffnm step6.0_minimization

# step6.1
gmx grompp -f step6.1_minimization.mdp -o step6.1_minimization.tpr -c step6.0_minimization.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx
gmx mdrun -deffnm step6.1_minimization
unset env GMX_MAXCONSTRWARN

# Equilibration
gmx grompp -f step6.2_equilibration.mdp -o step6.2_equilibration.tpr -c step6.1_minimization.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx
gmx mdrun -deffnm step6.2_equilibration
gmx grompp -f step6.3_equilibration.mdp -o step6.3_equilibration.tpr -c step6.2_equilibration.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx
gmx mdrun -deffnm step6.3_equilibration




# Copy back results from $SCRATCH to submission directory.
cd ..
cp -r gromacs $SLURM_SUBMIT_DIR/gromacs_out

exit 0

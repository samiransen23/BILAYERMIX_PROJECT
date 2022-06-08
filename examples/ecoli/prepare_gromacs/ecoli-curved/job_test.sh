mkdir -p out

# export paths of files that do not need changing
INPUT_PATH="/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-curved"

# modify topol.top
cp /Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-curved/topol.top .
echo '#include "toppar/martini_v2.0_lipids_all_201506.itp"' | cat - topol.top > temp && mv temp topol.top
echo '#include "toppar/martini_dummy.itp"'              | cat - topol.top > temp && mv temp topol.top

#run
# Minimization
## step6.0: softcore minimization
#gmx grompp -f mdp/step6.0_minimization.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c ${INPUT_PATH}/out_ecoli_centred.gro -o out/step6.0_minimization.tpr
#gmx mdrun -deffnm out/step6.0_minimization #-v
#
## step6.1: no softcore minimization
#gmx grompp -f mdp/step6.1_minimization.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c out/step6.0_minimization.gro -o out/step6.1_minimization.tpr
#gmx mdrun -deffnm out/step6.1_minimization #-v
#

## step 6.2: equilibration NVT timestep = 0.002 fs
#gmx grompp -f mdp/step6.2_equilibration.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c out/step6.1_minimization.gro -r out/step6.1_minimization.gro -o out/step6.2_equilibration.tpr
#gmx mdrun -deffnm out/step6.2_equilibration -rdd 1.3 -v


## step 6.3: equilibration NVT timestep = 0.005 fs
#gmx grompp -f mdp/step6.3_equilibration.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c out/step6.2_equilibration.gro -r out/step6.2_equilibration.gro -o out/step6.3_equilibration.tpr
#gmx mdrun -deffnm out/step6.3_equilibration -rdd 1.3 -v

# step 6.4: equilibration NVT timestep = 0.010 fs
#gmx grompp -f mdp/step6.4_equilibration.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c out/step6.3_equilibration.gro -r out/step6.3_equilibration.gro -o out/step6.4_equilibration.tpr
#gmx mdrun -deffnm out/step6.4_equilibration -rdd 1.3 -v

# step 6.5: equilibration NVT timestep = 0.015 fs
#gmx grompp -f mdp/step6.5_equilibration.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c out/step6.4_equilibration.gro -r out/step6.4_equilibration.gro -o out/step6.5_equilibration.tpr
#gmx mdrun -deffnm out/step6.5_equilibration -rdd 1.3 -v

# step 6.6: equilibration NVT timestep = 0.020 fs
#gmx grompp -f mdp/step6.6_equilibration.mdp -p topol.top -n ${INPUT_PATH}/index.ndx -c out/step6.5_equilibration.gro -r out/step6.5_equilibration.gro -o out/step6.6_equilibration.tpr
#gmx mdrun -deffnm out/step6.6_equilibration -rdd 1.3 -v

exit 0

#mkdir out

# Minimization
## step6.0: softcore minimization
gmx grompp -f ../mdp/step6.0_minimization.mdp -p topol.top -n index.ndx -c input_ionised.gro -r input_ionised.gro -o out/step6.0_minimization.tpr
gmx mdrun -deffnm out/step6.0_minimization -v -nsteps 400

# Equilibration
#gmx grompp -f ../mdp/step6.2_equilibration.mdp -p topol.top -n index.ndx -c out/step6.0_minimization.gro -r out/step6.0_minimization.gro -o out/step6.2_equilibration.tpr
#gmx mdrun -deffnm out/step6.2_equilibration -rdd 1.3 -nsteps 20000 -v
wait

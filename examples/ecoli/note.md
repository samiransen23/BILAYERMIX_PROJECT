Example run to generate curved double bilayer:
```bash
python BumPy/BUMPy/bumpy.py -f BumPy/BUMPy/examples/input_pdbs/small_flat_bilayer.pdb -s semicylinder_plane_U -z 10 \
-g r_cylinder:60 l_cylinder:100 r_junction:30 l_flat:200 l_margin:0  -o out.pdb \
-p topol.top -n index.ndx --gen_dummy_particles --dummy_grid_thickness 50
```

### E.Coli
### Method 1
#### Step 1: Build patch on CHARMM-GUI
A small flat bilayer patch was prepared on CHARMM-GUI with 67% POPE, 23% POPG and 10% CDL2 (cardiolipin) in a box of lengths `10.12598  10.12598   8.49522`.
It was neutralised with NA ions.  
Files are found in:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-flat/small/charmm-gui-5391855495`.  

#### Step 2: Equilibrate patch on Gromacs
Gromacs was run to energy-minimize (2 out of 2 minimisations that CHARMM-GUI provided) and equilibrate (2 out of the 5 equilbrations) according to the jobscript:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-flat/small/gromacs_out/job_charmmgui.sh`  
The final equilibrated flat bilayer patch is called: `step6.3_equilibration.*`

#### Step 3: Build required geometry with BUMPy
To build a U-shaped double bilayer:
Make a pdb from `step6.3_equilibration.gro` with selections `resname POPG or resname POPE or resname CDL2` and save as `input_flat_membraneonly.pdb`
The molecules are not necessarily whole, and BUMPy cannot handle broken molecules across pbc. So make them whole by:
```bash
gmx trjconv -f input_flat_membraneonly.pdb -pbc whole \
-s ~/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-flat/small/gromacs_out/step6.3_equilibration.tpr \
-n ~/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-flat/small/gromacs_out/index.ndx  \
-o input_whole.pdb
```

We can now use BUMPy to create different U-shaped geometries by suitably duplicating this patch.
For example:
```bash
python BumPy/BUMPy/bumpy.py -f input_whole.pdb -s semicylinder_plane_U -z 20 -g r_cylinder:60 l_cylinder:100 \
r_junction:0 l_flat:200 l_margin:0 -p topol.top -n index.ndx --gen_dummy_particles --dummy_grid_thickness 65 \
-o out_ecoli.pdb
```
Note: (1) `l_margin` is a new flag added personally to control a possible crack between the curved and flat part of the U-shape.  
&emsp;&emsp;&ensp;&ensp; (2) Optionally use  `--gen_dummy_particles --dummy_grid_thickness 50` flags to add dummy-particles. See documentation or publication for further info.  
&emsp;&emsp;&ensp;&ensp; (3) [BUMPy: A Model-Independent Tool for Constructing Lipid Bilayers of Varying Curvature and Composition](10.1021/acs.jctc.8b00765)

BUMPy generates box dimensions from the coordinates.
We want the structure to have extra space around it and also to be centred.  
We choose the output to be `gro`
```bash
gmx editconf -f out_ecoli.pdb -box 10 24 32 -o out_ecoli_centred.gro -c
```
#### Step 4: Run Gromacs
There is a bash script for it:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-curved/job_membrane_and_dummy.sh`  
Details:  
Minimize: 5000 steps  
NVT equilibrate: 10 - 50 ns broken down into different time-step steps.
NPT: Not meaningful because we are in vacuum.

<img src="/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-curved/images/semicylinder_plane_U/r_cylinder=60_l_cylinder=100.tga" width="500"/>


### Method 2
#### Step 1: Build patch on CHARMM-GUI
Same as Method 1
#### Step 2: Minimize patch on Gromacs
We only do minimisation and no NVT. Undesirable bumps were observed during NVT at this step in Method 1.
The final flat bilayer patch from this step is: `step6.1_minimization.*`
#### Step 3: Build required geometry with BUMPy
Same as method 1 but remember to replace `step6.3_equilibration.*` with `step6.1_minimization.*`.
Call BUMPy output `out_ecoli_U.pdb` and convert to `out_ecoli_U.gro`
#### Step 4: Complete the structure
1. We reflect the right U on the XY-plane. We also move the reflected structure slightly right (by 1.28nm) to match well the original.
We append it to the original:
```bash
gmx editconf -f out_ecoli_U.gro -o out_new.gro  -scale 1 1 -1 -translate 0 0 1.28
gmx editconf -f out_ecoli_U.gro -o 1.pdb
gmx editconf -f out_new.gro -o 2.pdb
cat 1.pdb 2.pdb > both.pdb
```
2. Some manual line deleting in `both.pdb`: Delete unnecessary lines like CRYST1 or END in the middle of the file.
3. Generate gro file, number residues correctly, center. You may remove the intermediate files no longer needed:
```bash
gmx editconf -f both.pdb  -box 10 24 64 -c -o both.gro
gmx genconf -f both.gro -o both_renumber.gro -renumber
rm out_new.gro 1.pdb 2.pdb both.pdb both.gro
mv both_renumber.gro out_ecoli_centred.gro
```
4. For Gromacs, we need an index file with `DUMMY` and `not_DUMY` particles in different groups:
```bash
gmx make_ndx -f out_ecoli_centred.gro -o index.ndx
```
>Select group as: `a 1 - 23895 | a 28285 - 52179`
```bash
sed 's/a\_1\-23895\_a\_28285\-52179/not\_DUMY/' index.ndx > index.ndx
```

5. For Gromacs, we also need an updated topology file:
We duplicate the content lines.

#### Step 5: Run Gromacs
There is a bash script for it:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-solvated-from10x10x8/job_membrane_and_dummy.sh`  
Details:  
Minimize: 5000 steps  
NVT equilibrate: 10 - 50 ns broken down into different time-step steps.
NPT: Not meaningful because we are in vacuum.

### Method 3
#### Step 1, 2, 3, 4
Same as method 2

#### Step 5: Solvate and ionise
`water.gro` is a box of Martini waters. `topol.top` is a copy of the topology of `out_ecoli_centred.gro`.  

Solvate:
```bash
method2buildinputpath="/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-complete-from-unrelaxed-patch"
cp ${method2buildinputpath}/topol.top .
gmx solvate -cp ${method2buildinputpath}/out_ecoli_centred.gro -cs water.gro 
-p topol.top -radius 0.21 -o ecoli_solvated.gro
```
We have specified a van der Waals radius of `0.21 nm` explicitly
because a deafult of 0.10 is too dense for Martini waters.
We also observed too much water within the lipid membranes and
these were removed with VMD with a suitable selection as:  
`not (resname W and within 6 of (name C1A or name C1A1 or name C1A2 or name C1B
or name C1B1 or name C1B2 or name C2A1 or name C2A2 or name C2B or name C2B1
or name C2B2 or name C3A or name C3B or name C4A or name C4A1 or name C4A2
or name C4B or name C4B1 or name C4B2 or name C5A1 or name C5A2 or name C5B1 or name C5B2)) `  
`ecoli_solvated.gro` was overwritten.
Also, modify the new (fewer) \# of `W` in `topol.top`

Ionize as:
```bash
gmx grompp -f ions.mdp -c ecoli_solvated.gro -p topol.top -o ecoli_solvated.tpr
gmx genion -s ecoli_solvated.tpr -o ecoli_ionised.gro -p topol.top -pname NA -neutral
gmx make_ndx -f ecoli_ionised.gro -o index.ndx
```

#### Step 6: Run Gromacs
Increase repulsive nature of the DUMY particles with C-tails
from `0.25810E-01` to `0.25810E+01`
You may need to use a rcouloumb, rvdw and rlist > 1.115.
I went ahead with the default.

There is a bash script for it:  `/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-solvated/job_membrane_and_dummy.sh`  
^ But I left this execution incomplete

### Method 4
#### Step 1: Build patch with 0% CA on CHARMM-GUI
Same as Method 5, except use APL=65 Angstrom on CHARMM-GUI

#### Step 2: Equilibrate patch on Gromacs
Attempt is to ensure no bumps on small patch, before going for the full structure.
The box relax to a lower APL=57.1 Angstrom (desired APL=61.5).  
Bumps almost don't appear!
The final equilibrated flat bilayer patch is called: `step6.4_equilibration.*`

#### Step 3: Build required geometry with BUMPy
Same as method 2 but remember to replace `step6.1_minimization.*` with `step6.4_equilibration.*`.

#### Step 4: Complete the structure
Same as method 2  
For `gmx make_ndx`:  
>Select group as:  
`a 1 - 24828 | a 29218 - 54045`
```bash
sed 's/a\_1\-24828\_a\_29218\-54045/not\_DUMY/' index.ndx > index.ndx
```
#### Step 5: Solvate and ionise

Same as method 3
Check the index file, and modify if needed.

#### Step 6: Run Gromacs
Same as method 3
Note that we run NPT now since we have water and ions.
Also, we have randomly stripped off some water, so NPT is mandatory.  
No bumps. Looks good, but it elongates forever.

Note that there is no CDL2 in this system, so I switched to a system with
2% CDL2 but starting from a larger box. Bumps appear again.
This means both an initial tight box as well as CDL2 are separately
responsible for the bumps.

We go to the next method (still no CDL2) but attempts to stop the forever
elongation.

### Method 5
#### Step 1-5
Same as Method 4

#### Step 6: Run Gromacs
##### Step 6.1  
To run NPT with dummy particles ---  
Increase DUMY position restraint to 10000.  
Use anisotropic Berendsen barostat for NPT until vacuum is removed.  
This initial equilibration ensures correct density of water
without changing the configuration of the lipids or dummy particles.

There is a bash script for it: `/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-from10.3x10.3x8.5/controlledelongation/step1_NPT`

##### Step 6.2  
To run NVT with dummy particles ---
Extract frame at 120ps where the vacuum is seen to have been removed as:
```bash
gmx trjconv -f gromacs_out/step6.2_equilibration.xtc -s gromacs_out/step6.2_equilibration.tpr -b 119 -e 120
-o step6.2_equilibration_120ps.gro
```
Restore dummy position restraint to 1000.  
Set mdp options to NVT.  
There is a bash script for it in: `/cluster/projects/nn4654k/samiran/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-from10.3x10.3x8.5/controlledelongation/step2_NVT`

##### Step 6.3a
To run NVT without dummy particles ---  
Extract last frame of previous step without dummy particles. Make corresponding
changes in topology and index files. These are automatically done by the jobscript.
Change `tc-grps` to `System`. Remove any mention of `DUMY` or `not_DUMY`.
There is a bash script for it in: `/cluster/projects/nn4654k/samiran/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-from10.3x10.3x8.5/controlledelongation/step3_NVT`

##### Step 6.4a
To run NPT without dummy particles ---
Use anisotropic Parinello-Rahman barostat.  
The system is stable. However some undesirable curvature is observed which 
may not have occured if there was an adhering surface, as is the goal.
So do an alternate step 6.3

##### Step 6.3b
To run NVT without dummy particles but with surface ---  
_Some building is needed first:_  
1. Extract last frame of previous step without dummy particles.
2. Make corresponding changes in topology and index files.  
3. Rotate frame so that xy-plane is the membrane-plane.
```bash
#steps 6.3b.1 - 3
STEP2_OUTPUT_PATH="/cluster/projects/nn4654k/samiran/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-from10.3x10.3x8.5/controlledelongation/step2_NVT/gromacs_out"
INPUT_PATH="/cluster/projects/nn4654k/samiran/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-complete-from10.3x10.3x8.5"
echo -e "13 \n q" | gmx trjconv -f ${STEP2_OUTPUT_PATH}/step6.6_equilibration.xtc -s ${STEP2_OUTPUT_PATH}/step6.2_equilibration.tpr -dump -1 -n ${INPUT_PATH}/index.ndx -o step6.6_equilibration_last.gro
sed '/DUMY.*/d' topol.top > topol_noDUMY.top
mv topol_noDUMY.top topol.top
echo -e "q" | gmx make_ndx -f step6.6_equilibration_last.gro -o index.ndx
gmx editconf -f step6.6_equilibration_last.gro -rotate 90 0 0 -box $(tail -1 step6.6_equilibration_last.gro | awk '{print $1 " "  $3 " " $2}') -translate 0 $(tail -1 step6.6_equilibration_last.gro | awk '{print $3}') 0 -o reoriented.gro
```
4. On VMD, make surface as `VMD > Extenstion > Modeling > NanotubeBuilder > Graphene Sheet`:  
Choose `(xedge, yedge) = Box length in after-rotation x,y` and `# layers = 3`. Call it `graphene.gro`.
5. Catenate `reoriented.gro` and `graphene.gro`. Do not center or scale (no `-c`, no `-scale`). 
Remember to renumber residues and change the box. Call it `both_renumber.gro`.
```bash
gmx editconf -f reoriented.gro -o 1.pdb
gmx editconf -f graphene.gro -o 2.pdb
cat 1.pdb 2.pdb > both.pdb
#manually delete extra lines in both.pdb
gmx editconf -f both.pdb -o both.gro
gmx genconf -f both.gro -o both_renumber.gro -renumber
rm 1.pdb 2.pdb both.pdb both.gro
```
6. Remove water inside the surface. In this case, the region was found to be `z<6.7`.
You can use VMD to save coordinates of selection. Call it `both_renumber_water_stripped.gro`
7. Add more water above than below (convention), to prevent the surface in the periodic image above from
pulling the lipid heads too.
8. Catenate `both_renumber_water_stripped.gro` and `extra_water_layer_translated.gro`.
Call it `input_step3b_NVT.gro`
```bash
#steps 6.3b.7 - 9
MARTINIWATER_PATH="/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-complete-solvated-from10x10x8"
gmx solvate -cs ${MARTINIWATER_PATH}/water.gro -box $(tail -1 step6.6_equilibration_last.gro | awk '{print $1 " " $2}') 6 -o extra_water_layer.gro
gmx editconf -f extra_water_layer.gro -translate 0 0 $(tail -1 step6.6_equilibration_last.gro | awk '{print $3}') -o extra_water_layer_translated.gro
gmx editconf -f both_renumber_water_stripped.gro -o 1.pdb
gmx editconf -f extra_water_layer_translated.gro  -o 2.pdb
cat 1.pdb 2.pdb > both.pdb
#manually delete extra lines in both.pdb 
lz=$(echo "$(tail -1 both_renumber_water_stripped.gro | awk '{print $3}') +6" | bc)
gmx editconf -f both.pdb -box $(tail -1 both_renumber_water_stripped.gro | awk '{print $1 " " $2}') ${lz} -o both.gro
gmx genconf -f both.gro -o both_renumber.gro -renumber
rm 1.pdb 2.pdb both.pdb both.gro translated.gro reoriented.gro \
graphene.gro both_renumber_water_stripped.gro \
extra_water_layer.gro extra_water_layer_translated.gro
mv both_renumber.gro input_step3b_NVT.gro
```
9. We decided (in an afterthought) to bring the surface closer, and remove water in the neighbourhood more wisely.
Use on `input_step3b_NVT.gro` a tcl script `surface_closer/surface_closer.tcl` which saves a new file `surface_closer/input_step3b_NVT.gro`.
```bash
source surface_closer/surface_closer.tcl
```
10. Change the box in Z by reducing it by about `0.5`.  
Rename the old `input_step3b_NVT.gro` as `surface_far.gro`.  
Renumber residues and modify index file.  
Modify topology.
```bash
#change box manually (coz fuck automating everything)
mv input_step3b_NVT.gro surface_far.gro
cd surface_closer
gmx genconf -f input_step3b_NVT.gro -o input_renumbered.gro -renumber
gmx make_ndx -f input_renumbered.gro -o index.ndx
#Select as follows: (change these)
#a 49657 - 99420
#a 173040 - 205456
sed 's/a\_49657\-99420/W/' index.ndx > temp.ndx
sed 's/a\_173040\-205456/W/' temp.ndx > temp2.ndx
mv temp2.ndx index.ndx
#modify topology manually (new water and ions)
```

Add neutralising ions (because we have removed some for the surface neighbourhood)
Again, renumber residues and modify index file.  
```bash
#add ions
gmx grompp \
-f ~/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-complete-solvated-from10x10x8/ions.mdp \
-c input_renumbered.gro -p topol.top -n index.ndx -o input_renumbered.tpr
gmx genion -s input_renumbered.tpr -o input_ionised.gro -p topol.top -pname NA -neutral -n index.ndx
#choose group 8 (containing 32417 atoms)

#modify index
echo -e "q" | gmx make_ndx -f input_ionised.gro -o index.ndx
```

_Some interactions with the surface need to be specified:_  
Surface particles GRA are chosen to be of their own type GRA.
POPG-head and W is of type P4,    
POPE-head and NA is of type Qd  
We choose GRA - P4 < 0, GRA - Qd > <0 and GRA - C1 > 0  
_We also use freeze groups for GRA in .mdp_

11. We observe crystallisation of water. As a second afterthought, we want the surface to be less dense.
We rebuild a system with such a surface in folder `surface_lessdense`
We change the lattice spacing to 3.2Ang from 1Ang.
We still use 3 layers. Save the surface as `surface.pdb`
Remove from `../input_ionised.gro` the surface and save as `nosurface.pdb`.
Combine the two to `both.pdb` and convert to gro:
```bash
gmx editconf -f both.pdb -o both.gro
```
Raise the surface by 1.5 using tcl and change the box in y-dimn to `63.30810` to ensure the full surface is inside
the box. Call it `input_step3b_NVT.gro` and update index:
```bash
gmx make_ndx -f input_step3b_NVT.gro -o index.ndx
rm 1.pdb 2.pdb both.pdb both.gro
```
Now we can run Gromacs. We have modified the `.mdp` files into NPT runs, since we need to get rid of the vacuum that
resulted from arbitrarily setting the box size.
The files are found in: `/cluster/projects/nn4654k/samiran/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete-from10.3x10.3x8.5/controlledelongation/step3b_NVT/surface_closer/surface_lessdense`

### Notes
> 1. For convenient VMD selection, a list of all C-like CG beads:
`name C1A or name C1A1 or name C1A2 or name C1B or name C1B1 or name C1B2 
or name C2A1 or name C2A2 or name C2B or name C2B1 or name C2B2 or name C3A 
or name C3B or name C4A or name C4A1 or name C4A2 or name C4B or name C4B1 
or name C4B2 or name C5A1 or name C5A2 or name C5B1 or name C5B2`


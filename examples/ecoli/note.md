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
#### Step 2: Equilibrate patch on Gromacs
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
>Select group as: `a 1 - 23895 | a 28285 - 52179`
```bash
gmx make_ndx -f both.gro -o index.ndx
sed 's/a\_1\-23895\_a\_28285\-52179/not\_DUMY/' index.ndx > index.ndx
```

5. For Gromacs, we also need an updated topology file:
We duplicate the content lines.

#### Step 5: Run Gromacs
There is a bash script for it:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-complete/job_membrane_and_dummy.sh`  
Details:  
Minimize: 5000 steps  
NVT equilibrate: 10 - 50 ns broken down into different time-step steps.
NPT: Not meaningful because we are in vacuum.

### Notes
> 1. For convenient VMD selection, a list of all C-like CG beads:
`name C1A or name C1A1 or name C1A2 or name C1B or name C1B1 or name C1B2 
or name C2A1 or name C2A2 or name C2B or name C2B1 or name C2B2 or name C3A 
or name C3B or name C4A or name C4A1 or name C4A2 or name C4B or name C4B1 
or name C4B2 or name C5A1 or name C5A2 or name C5B1 or name C5B2`

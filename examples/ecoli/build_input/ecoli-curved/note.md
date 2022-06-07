Example run to generate curved double bilayer:
```bash
python BumPy/BUMPy/bumpy.py -f BumPy/BUMPy/examples/input_pdbs/small_flat_bilayer.pdb -s semicylinder_plane_U -z 10 \
-g r_cylinder:60 l_cylinder:100 r_junction:30 l_flat:200 l_margin:0  -o out.pdb \
-p topol.top -n index.ndx --gen_dummy_particles --dummy_grid_thickness 50
```

### E.Coli
#### Build patch on CHARMM-GUI
A small flat bilayer patch was prepared on CHARMM-GUI with 67% POPE, 23% POPG and 10% CDL2 (cardiolipin) in a box of lengths `10.12598  10.12598   8.49522`.
It was neutralised with NA ions.  
Files are found in:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-flat/small/charmm-gui-5391855495`.  

#### Equilibrate patch on Gromacs
Gromacs was run to energy-minimize (2 out of 2 minimisations that CHARMM-GUI provided) and equilibrate (2 out of the 5 equilbrations) according to the jobscript:  
`/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/prepare_gromacs/ecoli-flat/small/gromacs_out/job_charmmgui.sh`  
The final equilibrated flat bilayer patch is called: `step6.3_equilibration.*`

#### Build required geometry with BUMPy
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
python BumPy/BUMPy/bumpy.py -f input_whole.pdb -s semicylinder_plane_U -z 20 -g r_cylinder:60 l_cylinder:100 r_junction:0 l_flat:200 l_margin:0 -o out_ecoli.pdb
```
Note: (1) `l_margin` is a new flag added personally to control a possible crack between the curved and flat part of the U-shape.  
&emsp;&emsp;&ensp;&ensp; (2) Optionally use  `--gen_dummy_particles --dummy_grid_thickness 50` flags to add dummy-particles. See documentation or publication for further info.  
&emsp;&emsp;&ensp;&ensp; (3) [BUMPy: A Model-Independent Tool for Constructing Lipid Bilayers of Varying Curvature and Composition](10.1021/acs.jctc.8b00765)

<img src="/Users/samiransen23/BILAYERMIX_PROJECT/examples/ecoli/build_input/ecoli-curved/images/semicylinder_plane_U/r_cylinder=60_l_cylinder=100.tga" width="500"/>

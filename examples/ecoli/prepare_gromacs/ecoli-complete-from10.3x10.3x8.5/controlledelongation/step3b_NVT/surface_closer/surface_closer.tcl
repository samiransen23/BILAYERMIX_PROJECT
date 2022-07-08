set gra [atomselect top "resname GRA"]
$gra moveby {0 0 8}
set strip [atomselect top "name W NA and ((z<50 and y<515 and y>108) or within 8 of resname GRA or z<19)"]
$strip moveby {0 0 1000}
set all [atomselect top "all"]
${all} moveby {0 0 -7.9}
set out [atomselect top "z<500"]
#$out writegro surface_closer/input_step3b_NVT.gro

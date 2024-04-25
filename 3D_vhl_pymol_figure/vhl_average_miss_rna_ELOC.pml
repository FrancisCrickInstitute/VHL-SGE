#For average function scores at each amino acid position for missense variant with RNA score higher or equal to -2: ---------------------------------
# Run in PyMOL prompt with the following command: run /Users/terwagc/Downloads/3D_vhl_pymol_figure/vhl_average_miss_rna_ELOC.pml

# get a white background
bg_color white

# Load the PDB file
load "/Users/terwagc/Downloads/3D_vhl_pymol_figure/1LM8_vhl_isolated.pdb", VHL
load "/Users/terwagc/PycharmProjects/vhl_dash_board/src/input/3d_structure/1lm8.cif", 1LM8


# Select and show only chain C
select C, 1LM8 and chain C
select elob, 1LM8 and chain B

deselect
color grey, all

# the desired view
set_view (\
     0.192830563,    0.356186360,   -0.914282560,\
     0.550915360,   -0.810347080,   -0.199504942,\
    -0.811971009,   -0.465226471,   -0.352485299,\
     0.000000000,    0.000000000, -227.297119141,\
    42.763023376,   19.123123169,  127.589080811,\
   156.116226196,  298.477966309,   20.000000000 )

zoom complete=1

# data2bfactor script obtained from here: http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/data2bfactor.py

# spectrumany script obtained from here: http://pymolwiki.org/index.php/Spectrumany
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/spectrumany.py

data2b_res VHL, /Users/terwagc/Downloads/3D_vhl_pymol_figure/Average_FS_at_AA_missense_only_rna_score_not_below_min2.txt
spectrum b, red_blue, VHL, minimum=-3.1, maximum=0.4



hide everything
show surface, VHL
show surface, C
show surface, elob

set transparency, 0, VHL
set cartoon_transparency, 0
set transparency, 0.75, C
set transparency, 0.75, elob

# Set color setting
set cartoon_transparency, 0
set depth_cue, 0  
set cartoon_oval_length, 1 
set cartoon_oval_width, 0.1 
set direct,0 
set ray_shadows, off
set orthoscopic, on
set light_count,5 
set spec_count,1 
set shininess, 3 
set specular, 0.2 
set reflect,1.3 
set ambient,0.3 
set backface_cull=0
set surface_quality, 1
set transparency_mode, 2

# for plastic shiny surface effect:
#set spec_reflect, 2.0
#set spec_power, 1500

set ray_opaque_background, on
ray 1000, 1000

# Save the image
png /Users/terwagc/Downloads/3D_vhl_pymol_figure/figures/Average_FS_at_AA_missense_only_rna_score_not_below_min2_eloc_TM1.png


# zoom resi 162+158 and VHL
set_view (\
     0.192830563,    0.356186360,   -0.914282560,\
     0.550915360,   -0.810347080,   -0.199504942,\
    -0.811971009,   -0.465226471,   -0.352485299,\
     0.000000000,    0.000000000,  -32.611091614,\
    52.278163910,   20.529214859,  115.635864258,\
  -23111.990234375, 23177.212890625,   20.000000000 )

hide everything
show surface, C
show surface, elob
show cartoon, C
show cartoon, elob
show cartoon, VHL
show stick, resi 158+162 and VHL


set ray_opaque_background, on
ray 1000, 1000

# Save the image
png /Users/terwagc/Downloads/3D_vhl_pymol_figure/figures/Average_FS_at_AA_missense_only_rna_score_not_below_min2_eloc_zoomTM1.png

## end ##
##For average function scores at each amino acid position for missense variant with RNA score higher or equal to -2: ---------------------------------
# Run in PyMOL prompt with the following command: run /Users/terwagc/Downloads/3D_vhl_pymol_figure/vhl_average_miss_rna.pml

# get a white background
bg_color white

# Load the PDB file
load "/Users/terwagc/Downloads/3D_vhl_pymol_figure/1LM8_vhl_isolated.pdb"
color grey, all

# the desired view
set_view (\
     0.746094227,   -0.518571973,   -0.417638302,\
    -0.120334201,    0.511894524,   -0.850579560,\
     0.654869497,    0.684867620,    0.319523156,\
     0.000000000,    0.000000000, -154.201324463,\
    42.209148407,   16.530614853,  132.441696167,\
  -942.812255859, 1251.214599609,  -20.000000000 )


# data2bfactor script obtained from here: http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/data2bfactor.py

# spectrumany script obtained from here: http://pymolwiki.org/index.php/Spectrumany
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/spectrumany.py

data2b_res 1LM8_vhl_isolated, /Users/terwagc/Downloads/3D_vhl_pymol_figure/Average_FS_at_AA_missense_only_rna_score_not_below_min2.txt
spectrum b, red_blue, 1LM8_vhl_isolated, minimum=-3.1, maximum=0.4

# make sure the figure is not truncated
zoom complete=1

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


set ray_opaque_background, on
ray 5400, 5400

# Save the image
png /Users/terwagc/Downloads/3D_vhl_pymol_figure/figures/Average_FS_at_AA_missense_only_rna_score_not_below_min2.png

## end ##
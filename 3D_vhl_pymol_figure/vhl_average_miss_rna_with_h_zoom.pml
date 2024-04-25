##For average function scores at each amino acid position for missense variant with RNA score higher or equal to -2: ---------------------------------
# Run in PyMOL prompt with the following command: run /Users/terwagc/Downloads/3D_vhl_pymol_figure/vhl_average_miss_rna_with_h_zoom.pml

# get white background
bg_color white

set_view (\
     0.650256097,   -0.727565527,   -0.218664601,\
     0.145568430,    0.401814997,   -0.904078186,\
     0.745635509,    0.556052804,    0.367194027,\
     0.000000000,    0.000000000,  -41.051048279,\
    34.559406281,   15.501226425,  142.261489868,\
  -33764.890625000, 33846.992187500,   20.000000000 )


# Load the PDB file
load "Downloads/3D_vhl_pymol_figure/1LM8_vhl_isolated.pdb", VHL
load "Downloads/3D_vhl_pymol_figure/1LM8_h_isolated.pdb", H


# data2bfactor script obtained from here: http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/data2bfactor.py

# spectrumany script obtained from here: http://pymolwiki.org/index.php/Spectrumany
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/spectrumany.py

data2b_res VHL, Downloads/3D_vhl_pymol_figure/Average_FS_at_AA_missense_only_rna_score_not_below_min2.txt
spectrum b, red_blue, VHL, minimum=-3.1, maximum=0.4

hide everything
show cartoon, VHL
#show surface, VHL
show sticks, H
color splitpea, H


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

# make sure the figure is not truncated
zoom complete=1

# display residues of interest as stick
show sticks, resi 111+115+117+88

# get high resolution image
set ray_opaque_background, on
ray 5400, 5400

# Save the image
png /Users/terwagc/Downloads/3D_vhl_pymol_figure/figures/Average_FS_at_AA_missense_only_rna_score_not_below_min2_HIF1A_zoom_a.png


# Zoom in on the pocket_resi selection with animation
zoom resi 111+115+117+88, animate=0
set depth_cue, 0  


# Save the image
set ray_opaque_background, on
ray 5400, 5400
png /Users/terwagc/Downloads/3D_vhl_pymol_figure/figures/Average_FS_at_AA_missense_only_rna_score_not_below_min2_HIF1A_zoom_b.png



## end ##
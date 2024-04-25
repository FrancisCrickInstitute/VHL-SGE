##For average function scores at each amino acid position for missense variant with RNA score higher or equal to -2: ---------------------------------
# Run in PyMOL prompt with the following command: run /Users/terwagc/Downloads/3D_vhl_pymol_figure/red_blue_spectrum_bar.pml

# get a white background
bg_color white

# spectrumbar script obtained from here: https://pymolwiki.org/index.php/Spectrumbar
run /Users/terwagc/Downloads/3D_vhl_pymol_figure/spectrumbar.py

# Get spectrum bar -----------------------------------
spectrumbar red, blue

# make sure the figure is not truncated
zoom complete=1

# Set color setting
set depth_cue, 0  
set cartoon_oval_length, 1 
set cartoon_oval_width, 0.1 
set ambient,0.3 
set ray_shadows, off
set orthoscopic, on
set direct=0.7
set reflect=0.0
set backface_cull=0
set surface_quality, 2
set transparency_mode, 2

# for plastic shiny effect:
#set spec_reflect, 2.0
#set spec_power, 1500


# Set the white background
set ray_opaque_background, on
ray 5400, 5400
png /Users/terwagc/Downloads/3D_vhl_pymol_figure/figures/red_blue_spectrum_bar.png

## end ##
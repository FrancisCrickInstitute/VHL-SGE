All isolated pdb files used are generated using the isolate_chains.pml file (from the VHL dash board repo)

To generate the VHL surface with average missense variants score for variant with RNA_score >= -2 :
1. Open the vhl_average_miss_rna.pml with a text editor, modify all paths with your local path and save the file.
2. Open the PyMol application 
3. In the PyMol prompt, run the following command:
PyMOL> run <your_path>/3D_vhl_pymol_figure/vhl_average_miss_rna.pml

To add the HIF 1A, repeat the steps 1 to 3 with the file vhl_average_miss_rna_with_h.pml
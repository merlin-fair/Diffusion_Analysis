# Diffusion_Analysis
Scripts for pre-processing and analysis of diffusion MRI data. 
This script contiain the stepts for obtaining the tractography of **n** sujects of your interest, and it also proportionate the .csv document that helps to construct a connectivity matrix / connectome.  
To see the matrix it's necessaire to change the programattion language, from linux/bash to MATLAB. Here it's shown too the code for visualizate your connectome in Matlab. 

## How does this work? 
Here it's shown a process that works with two scripts: a script-parent and a script-child. Both work with bash.
In the script-parent, you'll have first to write the absolute direction (parent) where your subjects are. Then, you have to place manually the exact names of the subjects that you want to process. And finally, you can write the rute where your script-child is. 
* You'll have to be careful because there are some process that take their time, so maybe it's a good idea to make a test with few sujects at first (it depends on your computer's / cluster's power). 

This script uses the enviroments of MRtrix, FSL, ANTs and Freesurfer. At some point inside the pipeline, you can find a part of Python inside the bash script. At the end of this script, 
and also it uses the conda enviroment. 

Requirements before running the code: 
Establish the following modules on the linux terminal. Here are written the specifiq versions that make this code work but it may change at some point:
* fsl/6.0.7.22
* ANTs/2.4.4
* MRtrix/3.0.4
* freesurfer/7.4.1.
  
Install the following conda / Python libraries: 
* 

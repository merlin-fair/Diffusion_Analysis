# Diffusion_Analysis
Scripts for pre-processing and analysis of diffusion MRI data. 
This script contiain the stepts for obtaining the tractography of __n_ sujects of your interest, and it also proportionate the .csv document that helps to construct a connectivity matrix / connectome. On the MATLAB part, you can use this code to make a simple group analysis. 

## How does this work? 
Here it's shown a pipeline that works with two scripts: a script-parent and a script-child. Both work with bash.

In the script-parent, you'll have first to write the absolute path (variable parent) where your subjects are. Then, you have to place manually the exact names of the subjects that you want to process. And finally, you can write the rute where your script-child is. 
* You'll have to be careful because there are some process that take their time, so maybe it's a good idea to make a test with few sujects at first (it depends on your computer's / cluster's power) and see how it goes.

The child script contains all the commands that will correct the distorssions of you images, and at the end of the process, you'll get the tractography that you're looking for. Also, it gives you a .csv document which you can use to obtain a connectivity matrix.

If you want to visualize the matrix, you can use the MATLAB script. You will be able to see the matrices / connectomes of all your sujects if you want, but also you can access to a script that helps you to get a simple analysis between groups. 

This script uses the enviroments of MRtrix, FSL, ANTs and Freesurfer. At some point inside the pipeline, you can find a part of Python inside the bash script. 

## How do I organise my data? 
This pipline was thought as a tool to get the analysis between groups. For that reason, two types of data was used: the one for the control-group and the other for the experimental-group. 
* You'll need to have a folder in the linux terminal (main folder). Inside of this, you'll have a folder for the control-group and other for the experimental-group.
* You should move your corresponding subjects on each folder.
* Each subject should contain the next files to make the script work: 
*   dwi.bval
*   dwi.bvec
*   dwi.jason
*   dwi.nii.gz
*   T1w.json
*   T1w.nii.gz
*   Just an advice, you can make an __original_files_ folder inside of each subject. The original_files will contains the same files just mentioned. This is a good practice because you can access in an easy way to the principal files of your subjects after you get all the new files from the tractography. Something like this:
*   Subject 1
*   ______ All the
*   ______ files just
*   ______ mentioned
*   ___________________  Folder with the __original_files_

If you want to continue with the MATLAB pipeline, you should organise your data like this: 
* Once you get the results of all your subjects, you have to take the file named __connectome_sift.csv_ of each subject.
* Inside the MATLAB terminal, you should make a very similar organisation as in Linux: you'll have your main folder, and inside, a control-group folder and an experimental-group folder. You have to paste the file __connectome_sift.csv_ of each subject as corresponding.
* Also, there is adocument named __Copy_of_label_fs_default.txt_ which contains the information to write the corresponding labels of the connectome. You should colocate this on your main folder. 

## Requirements before running the code: 
Establish the following modules on the linux terminal. Here are written the specifiq versions that make this code work but it may change at some point:
* fsl/6.0.7.22
* ANTs/2.4.4
* MRtrix/3.0.4
* freesurfer/7.4.1.
  
Install the following conda / Python libraries: 
* import sys
* import os
* import pandas as pd
* import numpy as np

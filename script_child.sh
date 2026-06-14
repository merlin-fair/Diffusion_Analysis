#!/bin/bash
#Este es el script hijo, que permite llevar a cabo las tractografías. Obecede al script padre, que le indica qué sujetos procesará.
#Para copiar los datos de una carpeta en sí misma hacer: cp -r sub-10217/ sub-10217/original_files

set -x

------------------------------------------------------------------------------
# 1. CONFIGURACIÓN DE RECURSOS
# ------------------------------------------------------------------------------

#Aqui el primer argumento $1 que viene del script padre se guarda en la variable subject_dir.  
subject_dir="$1"
N_THREADS=6  # Como Opuntia tiene 24 núcleos NCOR, 6 hilos permite correr 3 sujetos a la vez.

# -z significa una lista vacía. Entonces, si subject_dir no tiene ningun argumento, en este caso, ninguna carpeta, entonces acaba el proceso. 
if [ -z "$subject_dir" ]; then 
	echo "Error: no se proporcionó el directorio de sujeto"
	exit 1 
fi 

echo "Procesando: $subject_dir"
	
cd "$subject_dir" || { echo "No se pudo entrar a $subject_dir"; exit 1; }

# ------------------------------------------------------------------------------
# 2. CONTROL DE CALIDAD PREVENTIVO
# ------------------------------------------------------------------------------
# Comprobamos archivos necesarios antes de empezar el proceso de 3 horas

if [[ ! -f "dwi.nii.gz" || ! -f "dwi.bvec" || ! -f "dwi.bval" ]]; then
    echo "ERROR: Faltan archivos dwi/bvec/bval en $subject_dir. Abortando."
    exit 1
fi

#------------------------------------------------------------------------------
# 3. PIPELINE DE PROCESAMIENTO (MRtrix3 + FSL)
# -----------------------------------------------------------------------------
echo "-----------Comienza proceso de denoising---------"
echo "---Conversión a .mif---"
mrconvert dwi.nii.gz dwi.mif -fslgrad dwi.bvec dwi.bval -nthreads $N_THREADS -force
echo "---Aplicar el denoising a la imagen---"
dwidenoise dwi.mif dwi_den.mif  -noise noise.mif -nthreads $N_THREADS -force
echo "--- proceso tardado: dwifslpreproc (eddy) ---"
dwifslpreproc dwi_den.mif dwi_den_preproc.mif -nocleanup -pe_dir AP -rpe_none -nthreads $N_THREADS  -eddy_options " --slm=linear --data_is_shelled" 
dwibiascorrect ants dwi_den_preproc.mif dwi_den_unbiased.mif -bias bias.mif -nthreads $N_THREADS -force
dwi2mask dwi_den_unbiased.mif mask.mif -nthreads $N_THREADS -force

echo "--- Iniciando pipeline de tractografia ---" 

#Cambiar el nombre de la T1 
mv T1w.nii.gz t1.nii.gz

if [[ ! -f "t1.nii.gz" ]]; then
    echo "ERROR: Faltan el archivo t1.nii.gz en $subject_dir. Abortando."
    exit 1
fi 

dwi2response dhollander dwi_den_unbiased.mif wm.txt gm.txt csf.txt -voxels voxels.mif 
	#Para empezar a crear ls FODS 
	dwi2fod msmt_csd dwi_den_unbiased.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif 
	mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif 
	#Llevar a cabo normalizacion para analisis futuros 
	mtnormalise wmfod.mif wmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif
	

echo "=== AQUI VA A EMPEZAR PYTHON ==="

# ==============================================================================
# 1. CONFIGURACIÓN DEL ENTORNO (¡Aquí van las líneas!)
# ==============================================================================
export FSLDIR=/home/inb/soporte/lanirem_software/fsl_6.0.7.22
export PATH=${FSLDIR}/bin:${PATH}
. ${FSLDIR}/etc/fslconf/fsl.sh

echo "=== Ejecutando Python dentro de Bash ==="

/home/inb/lailsongm/miniconda3/envs/servicio/bin/python << 'EOF'

from dipy.io.image import load_nifti, save_nifti
import numpy as np
import sys
sys.path.append('/misc/opuntia/lailsongmiranda')
import preprocessing as ppr
from subprocess import run
from os.path import join
from os import listdir as ls

#Empieza el código
run(['fslreorient2std', 't1.nii.gz', 't1_std.nii.gz'])
run(['first_flirt', 't1_std.nii.gz', 't1_std_mni.nii.gz'])

mathexa = np.loadtxt('t1_std_mni.mat')
np.savetxt('t1_std_mni_dec.mat', mathexa, delimiter=' ', fmt='% f')

caudadoleft = '/home/inb/soporte/lanirem_software/fsl_6.0.7.4/data/first/models_336_bin/L_Caud_bin.bmv'
caudadoright = '/home/inb/soporte/lanirem_software/fsl_6.0.7.4/data/first/models_336_bin/R_Caud_bin.bmv'
dir_modelos = '/home/inb/soporte/lanirem_software/fsl_6.0.7.4/data/first/models_336_bin/05mm/'
model05 = ls(dir_modelos)
model05

model05.remove('R_Cereb_05mm.bmv')

model05.remove('L_Cereb_05mm.bmv')

modelpaths = [join(dir_modelos, x) for x in model05]

modelpaths.append(caudadoleft)

modelpaths.append(caudadoright) 

modelpaths

modelosnii = []
for y in modelpaths: 
    nombre = y.split('/')
    nombre_bmv = nombre[-1]
    nombrenii = nombre_bmv.split('.')[0]+'.nii.gz'
    modelosnii.append(nombrenii)
    run(['run_first', '-i', 't1_std.nii.gz', '-v', '-t', 't1_std_mni_dec.mat', '-n', '40', '-o', nombrenii, '-m', y ])

listadata = []
mdata, mafin = load_nifti(modelosnii[0])
shapedata = mdata.shape
shapedata

submask = np.zeros(shapedata)
np.any(submask)

for strmodel in modelosnii:
    datastruct, afinstruct = load_nifti(strmodel)
    submask = np.where(datastruct, 1, submask)
    
#submask = np.zeros(shapedata)
#np.any(submask)
save_nifti('subcort_mask.nii.gz', submask.astype(np.float32), afinstruct)

run(['bet', 't1_std.nii.gz', 't1_std_bet', '-f', '0.25', '-m', '-B'])

run(['fast', 't1_std_bet.nii.gz'])

csfdata, csfafin = load_nifti('t1_std_bet_pve_0.nii.gz') 
gmdata, gmafin = load_nifti('t1_std_bet_pve_1.nii.gz') 
wmdata, wmafin = load_nifti('t1_std_bet_pve_2.nii.gz') 
threet = [csfdata, gmdata, wmdata]

whtsubcort = []
for tissue in threet:
    nosub = np.where(submask, 0, tissue)
    whtsubcort.append(nosub)
    
csf, gm, wm = whtsubcort

patolog = np.zeros(shapedata)
fivetiss = np.array([gm.T, submask.T, wm.T, csf.T, patolog.T])
save_nifti('fivearray.nii.gz', fivetiss.T, csfafin)    

EOF

echo "=== CONTINUAMOS CON BASH ==="

mrconvert fivearray.nii.gz fivearray.mif 
	dwiextract  dwi_den_unbiased.mif  - -bzero | mrmath - mean mean_b0.mif -axis 3 
	fslroi fivearray.nii.gz 5tt_vol0.nii.gz 0 1
	mrconvert mean_b0.mif mean_b0.nii.gz
	flirt -in mean_b0.nii.gz -ref 5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat diff2struct_fsl.mat
	transformconvert diff2struct_fsl.mat mean_b0.nii.gz fivearray.nii.gz flirt_import diff2struct_mrtrix.txt
	mrtransform fivearray.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif 
	5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif

	#PARA EMPEZAR A CREAR LAS STREAMLINES 
	tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmSeed_coreg.mif -nthreads 8 -maxlength 250 -cutoff 0.06 -select 10000000 wmfod_norm.mif tracks_10M.tck
	#Si queremos visualizar el output tracks_10M.tck 
	#tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck
	#mrview den_unbiased.mif -tractography.load smallerTracks_200k.tck
	tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt -nthreads 8 tracks_10M.tck wmfod_norm.mif sift_1M.txt
  
# ==============================================================================
# 5. COMANDO PESADO RECON-ALL 
# ==============================================================================

#Este es el código para correr RECON-ALL 
echo "Comenzando proceso de recon-all"

name=$(basename "$subject_dir") #extrae el nombre de la carpeta y lo guarda en la variable nombre 
subject_name="reconal_${name}"

export SUBJECTS_DIR="$subject_dir" 
	pwd 
	mkdir -p "$SUBJECTS_DIR"  #aqui se crea la carpeta de salida. -p significa 'crea la carpeta si no existe'. si existe no da error. 
	
if [ ! -f t1.nii.gz ]; then 
	echo "No existe t1.nii.gz en $subject_dir"
	exit 1 
fi

echo "Lanzando recon-all para $subject_name" 
recon-all -s "$subject_name" -i t1.nii.gz  -all

OUT_DIR="$subject_dir/connectome"

echo " subject dir es "$subject_dir" "

mkdir -p "${OUT_DIR}"
chmod +x "${OUT_DIR}"

labelconvert "$subject_dir"/"reconal_${name}"/mri/aparc+aseg.mgz "$FREESURFER_HOME"FreeSurferColorLUT.txt /home/inb/soporte/lanirem_software/mrtrix_3.0.4/share/mrtrix3/labelconvert/fs_default.txt "$OUT_DIR"/nodes_t1.mif -force

mrtransform "${OUT_DIR}"/nodes_t1.mif -interp nearest -linear "$subject_dir"/diff2struct_mrtrix.txt -template "$subject_dir"/mean_b0.mif "${OUT_DIR}"/nodes_dwi.mif -force

tck2connectome "$subject_dir"/tracks_10M.tck "${OUT_DIR}"/nodes_dwi.mif "$subject_dir"/connectome_sift2_"${name}".csv -tck_weights_in "$subject_dir"/sift_1M.txt -out_assignments ${OUT_DIR}/assignments.txt -symmetric -zero_diagonal -assignment_radial_search 2 -force


   

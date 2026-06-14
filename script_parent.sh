#!/bin/bash
#Este es el script padre que tendrá las rutas absolutas y relativas de los sujetos y del script hijo.

set -x
#En parent, se debe de copiar la dirección de la carpeta donde están los sujetos controles o experimentales.
parent="/misc/your_computer/main_folder_control"

#En la primera línea del for irán los nombres iguales de los sujetos que se quieren someter al proceso. Estos sujetos deben de estar dentro de la ruta parent para que funcione el script.
for suj in sub-1 sub-2 sub-3 do;  
	target_dir="$parent/$suj"
	if  [ -d "$target_dir" ]; then 
	echo "Lanzando proceso para: $suj"
	fsl_sub -N test_${suj} bash /misc/your_computer/main_folder/script_child.sh "$target_dir"
	else 
	echo "Error: no se encontro el directorio $target_dir"
	fi
done

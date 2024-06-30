#!/bin/bash

#######################################
# python script
script_annotation="Gene_Hits_Annotation_Args.py"
#Folder uniref IDs and kegg are 
base_exec_folder="/mnt/f/Cesar_Tesis/Hits_To_UniRef_KOs" # Se debe ejecutar desde aqui ya que es donde estan los unirefIDs
################Folders of all RUNs
folde_input='' # Place where your files are
diamond_Rxs='/mnt/f/Cesar_Tesis/MGWAS_Textile/Metagenoma_Completo/MapeoENS_Bowtie2/Databases/PositiveMG_DB/DB2/Diamond_Annotation/DB2_allproteins.tsv'
folder_output_name= # Change accoring to your folder
####################################### 
# Process
for file in Significan_Gene_Hits_1.txt Significan_Gene_Hits_2.txt Significan_Gene_Hits_3.txt; do # Change according to your files
    # Folder de ejecucion
    cd "$base_exec_folder" || exit 1
    python $script_annotation -i "$folder/$file" -b "$diamond_Rxs" -k "KOs.keg" -o "$folder/$folder_output_name/$file"
done

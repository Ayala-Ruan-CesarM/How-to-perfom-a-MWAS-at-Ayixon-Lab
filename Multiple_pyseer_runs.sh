#!/bin/bash
# Fill the values with the file path
phenotypes=''
unitigs=''
lmm_model=''
covariates=''

#Scripts
script_qplot='/mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/Scripts_Analisis_Datos_Pyseer/QQplot_only.py'

# Initial values for the cycle  
min_af=0.05
OUTPUT_FOLDER_BASE="Run"
# Loop for the desired number of iterations
while (( $(echo "$min_af <= 0.30" | bc -l) )); do
    # Construct the run14 variable
    OUTPUT="${OUTPUT_FOLDER_BASE}${min_af}"
    
    # Ensure the output directory exists
    mkdir -p "$OUTPUT"

    # Run the command
    pyseer --lmm --phenotypes "$phenotypes" \
           --kmers "$unitigs" \
           --load-lmm "$lmm_model" \
           --covariates $covariates \
           --use-covariates 5q 6q 7q 8q 12q 13q \
           --min-af "$min_af" \
           --print-samples \
           --output-patterns "$OUTPUT/Kmer_patterns.txt" \
           --cpu 40 > "$OUTPUT/Association_Result_kmers.txt"
    
    # Quantile-Quantile plot 
    cd "$OUTPUT/"
    if [ -f "Association_Result_kmers.txt" ]; then
        python $script_qplot -i Association_Result_kmers.txt
    else
        echo "Association result file does not exist. Re run the model"
    fi
    cd ..
    # Increase the --min-af value by 0.05
    min_af=$(echo "$min_af + 0.05" | bc)
done

# How-to-perfom-a-MWAS-at-Ayixon-Lab
This repository intens to serve as a blueprint to perform a microbiome-wide association study on the archezoa server.
The main software useded is Pyseer (Lees et al., 2017), how ever others are also used, if not mentioned all are is installed on $PATH.
For further Pyseer information refer to : https://pyseer.readthedocs.io/en/master/tutorial.html 

## Prepare input files

We'll assume that all assemblies are on a folder, hereinafter referred to as "working_folder".
* TIP:
When working with several files and if they have a binary clasificaition e.g. Positive/Negative, I recommed to rename each sample according to its classification, that will help to keep a good track.
Additionally, add to the name a unique pattern or a numeric series for example, either all the positive samples has the word "Pos" on its name or all positive samples goes from 1 to X.
The same should be done for the negative samples.

Create a list that segregates the names of each group sample.
```
cd $working_folder/
ls *Pos* > Positive_samples.list
ls *Neg* > Negative_samples.list
```

First, a "phenotype" file is requiere. This is a tab separeted file with two columns, header requiere.
When working with a continuos phenotypes the Phenotype.pheno file is still requiere but the second column is now numerical
```
# Add binary outcome to samples
awk '{print $0"\t"1}' Positive_samples.list > Positive.temp
awk '{print $0"\t"0}' Negative_samples.list > Negative.temp

# Concatenate and add header (Sample and Phenotype) customizable. 
echo -e "Sample\tPhenotype" | cat - Positive.temp Negative.temp > Phenotype.pheno

# Remove the extension from the .pheno file
sed s'/.fasta//'g

# Clean up
rm -f *.temp
```

Second, a variant file calling is requiere. This file is used to optain the unitigs
If no RAM constraints are met:
```
ls -d -1 $PWD/*.fasta > Unitigs_input.txt
```
¿How to know that apriori? 
If your dataset has more than 40 GB of information... don't even bother. 
With a dataset of 20-30 Gb you are good to go. 
Otherwise, the solution I found to this problem was to generate the unitigs in batches for that you'll need:
A file that has the absolute path of each group sample or a evenly distributed files. 
```
sed "s|^|$PWD/|" Positive_samples.list > Positive_samples_unitigs_input.txt
sed "s|^|$PWD/|" Negative_samples.list > Negative_sampless_unitigs_input.txt
```
## Variant (unitigs) Calling
We are going to assume the RAM constraints scenario, for further information visit [unitig-caller](https://github.com/bacpop/unitig-caller)  
Therefore we are going to use the --call and --query options as follows:
```
unitig-caller --call --refs Positive_samples_unitigs_input.txt --pyseer --threads 4 --out Positive_unitigs
```
The output file contains all the unitigs present only in the samples indicated by the --refs option.  
The .pyseer extension is added by default.  
Now, let's querry those unitigs onto the remaining samples:  
```
unitig-caller --query --ref Negative_sampless_unitigs_input.txt --unitigs Positive_unitigs.pyseer --pyseer --threads 4 --out Positive_Negative_unitigs
```
This time the output file contains the unitigs from the Positive samples that are found in the negative samples.  
However, as unitig-caller was not designed to run in batches, a few lines in the file are corrupt and needs to be removed.  
This bug was catched while running Pyseer as it prints to screen unitigs that
The correct format of each unitig line must be:  
> ATGACATGACATGACATGACATGACATGACT | Sample_1:1 Sample_2:1 Sample_3:1 Sample_X:1

It is highly advise to create a back-up file and to keep a count of the number of lines in both the original and corrected file.  
```
cp Positive_Negative_unitigs.pyseer Positive_Negative_unitigs_backup.pyseer ; gzip -v -9 Positive_Negative_unitigs_backup.pyseer
```
To fix the file in one line:  
```
sed -i -e 's/|$//' -e 's/| //2' Positive_Negative_unitigs.pyseer ; grep -v -E '\|\s*$' Positive_Negative_unitigs.pyseer > temp_file && mv temp_file Positive_Negative_unitigs.pyseer
```
## Desining a Kinship Matrix 
I explored three ways to include the relashionship bewteen metagenomic samples in the LMM those are:  
* Phylogenetic:  
Compute a phylogenetig tree from your samples.  
```
bash /mnt/f/Cesar_Tesis/MGWAS_Textile/JolyTree2.0/JolyTree/JolyTree.sh -i $working_dir -b MetagenomeTree -t 30
```
The "-s" is set automatically as it is calculated as the  proportion (up to 1.00) of the average genome size  
but in my observatiosn the higher the better.  
To keep a track of the current sketch size the output print to the screen can be written to a file like this:  
```
{ /usr/bin/time -v /mnt/f/Cesar_Tesis/MGWAS_Textile/JolyTree2.0/JolyTree/JolyTree.sh -i $working_dir -b MetagenomeTree -t 30 ; } > Screen_Output.txt 2>&1
```
```
grep "sketch size:" Screen_Output.txt 2>&1
```
The output number is the actual sketch size used by JolyTree which can be increase.  
Now, the MetagenomeTree.nwk output from Jolytree is the input to create the phylogenetic matrix with Pyseer script.  
```
python /mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/pyseer_scripts/phylogeny_distance.py --lmm MetagenomeTree.nwk > Phylogeny_similarity_matrix.tsv
```
```
# Remove the extension from the similarity matrix file
sed -i s'/.fasta//'g Phylogeny_similarity_matrix.tsv
```
* Genotype matrix:  
For this will need calculate a design matrix of variant presence absence to calculate the kinship matrix.  
Which can be calcualted from the unitigs file or a gene presence absence (rtab) or vcf file. However, high attención is requiere in order
to avoid "dilution" phenomena. That is exclute the genetic variants used to create Genotype matrix from the association model.  
We'll use a python script form Pyseer. 
```
awk '{print $1}' Phenotype.pheno > sample_list.txt
```
```
similarity_pyseer --kmers Positive_Negative_unitigs.pyseer sample_list.txt > Genotype_kinship.tsv
```
* Distance Metric:    
You will need a ecological beta diversity metric with size NxN and the R package "MiRKAT".  
The script Generar_KernelMatrix_fromBray.R from this project [auxiliary scripts](https://github.com/Ayala-Ruan-CesarM/Dye_MWAS_Aux_Scripts) is need it.  
The execution is stated in that repository.  
From my project this strategy was the most promising. However, that doesn't imply that's the best for your data.
## Performing and linear mixed model on Pysser.
So far up to this point we create the three inputs requiere for Pyseer.  
**Phenotype file**: Phenotype.pheno  
**Variant file**: Positive_Negative_unitigs.pyseer  
**Kinship matrix**: Similarity_matrix.tsv  
Now, we are set to perform a association test on the data as follwos:
```
pyseer --lmm --phenotypes Phenotype.pheno --kmers Positive_Negative_unitigs.pyseer --uncompressed \
--save-lmm model.npz --similarity Similarity_matrix.tsv --covariates Covariates_file.tsv \
--use-covariates 2 3 6q 7q 8q --print-samples --output-patterns Kmer_patterns.txt --cpu 30 > Association_Result_kmers.txt
```
* TIP: The "--save-lmm" option is used to save to a .npz file the deconstruction of the similarity matrix.
In subsequent runs the model is load with the option "--loadl-lmm" and the "--similarity" option is no longer required.

In order to include covariates on to the associaiton model a tab separeted file is required. The first column contains the Sample name as in the phenotype file.  
The other columns (with header) contains the information. Pyseer accepts both quantitative and qualitative covariates. 
The "--use-covariates" option tells Pyseer whih covariate use according to this column position. If the covariate is quantitative the "q" letter is requiered.  
The "--print-samples" adds to the "Association_Result_kmers.txt" file a column indicating in which sample each variant is and isn't.  

To perform several Pyseer Runs increasing the "--min-af" in 0.05 each time, use the Multiple_pyseer_runs.sh script provided.  
Copy the script to your working directory and add the path to the input files then execute: 
```
bash Multiple_pyseer_runs.sh
```

## Variant interpretation and annotation.
This secction describes a basic interpretation of the "best" pyseer results according to what the previous QQ-plots tells and what I used in my project.  
Additional information is in [Pyseer_tutorial.](https://pyseer.readthedocs.io/en/master/tutorial.html).  
Here, I am using the full path in the archezoa server but the scripts can be also be found in [auxiliary scripts](https://github.com/Ayala-Ruan-CesarM/Dye_MWAS_Aux_Scripts).  

First, get the significance threshold.  
```
python /mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/pyseer_scripts/count_patterns.py Kmer_patterns.txt
```
Which will output something like:
> Patterns:       331026
> Threshold:      1.51E-07

Second, assest the distribution of a the variant effect size in the result with a Volcano Plot.  
Each point in the graph correspons to a variant and it's possition is located in the **X-axis** according to the effect size (Beta), in the **Y-axis** according to the LRT-pvalue
and the **color** is the standard error of the effect size. This also help to assess the model fitting. 
```
python /mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/Scripts_Analisis_Datos_Pyseer/Volcano_plot.py -i Association_Result_kmers.txt -o Figure_name -t 1.51E-07
```
Third, as we are interested only in the genetic variants that actually are associated with the phenotype (Effect size > 0.0).  
From the Association_Result_kmers.txt file, let's extract those variants and save it to another file.  
```
awk '$5>0.0 {print $0}' Association_Result_kmers.txt > Variants_Beta_Positive.txt
```
Fourth, unlike the official recommendation from pyseer documentation, we do the annotation over all the variants and latter over the gene hits the threshold is applied.  
```
python /mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/pyseer_scripts/kmer_mapping/annotate_hits.py Variants_Beta_Positive.txt DB_Reference.txt Annotated_Variants_Veta_Positive.txt
```
The "annotate_hits.py" was modified to be able to use 30 threads.  
The DB_Reference.txt is the same as describe in the original documentation.
* TIP: If you want to use the Textil Dye database as reference use the file place at: /mnt/f/Cesar_Tesis/MGWAS_Textile/Metagenoma_Completo/MapeoENS_Bowtie2/Databases/ReferenceD2.txt


Fifth, generation of the gene hits file. And summarizaiton.    
```
python /mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/pyseer_scripts/summarise_annotations.py Annotated_Variants_Veta_Positive.txt > Gene_hits.txt
```
```
python /mnt/f/Cesar_Tesis/MGWAS_Textile/scripts/Scripts_Analisis_Datos_Pyseer/Script_to_sumarize.py Gene_hits.txt
```
What "Script_to_sumarize.py" script does is to look whether a gene ID is repeated, if all but the "hits" column has the same values it merge both annotations and adds up the total number
of hits. For example if there are two "rpoC_12" genes with 20 and 2 hits and all other statistics are the same, the rpoC_12 updates to 22 hits.  
Finally, the scripts keeps the original file, the output keeps the name but add a "Sum" at the beggining of the file name "SumGene_hits.txt".  

Sixth, filter by LRT-pvalue and generate the significant gene hits file:  
```
awk '$3>1.51E-07 {print $0}' SumGene_hits.txt > Significant_Gene_Hits.txt
```

## Annotationg hypothetical proteins with Uniref90 database on archezoa server
From the Significant_Gene_Hits.txt we are going to extract the hypothetical proteins according to it's common pattern ID. This is dependent of the database used for the variant annotation.
Let's assume that the initial annotation of the database or geneme used as reference was made with PROKKA. So we are going to assume that all the hypotethical proteins has the common pattern "HYPO" 
and that the .faa (amino acid) file exist.  
```
awk '{print $1}' Significant_Gene_Hits.txt | grep "HYPO" > Protein_Ids.txt
```
Now, from the .faa file.  
```
seqkit grep -f Protein_Ids.txt -o Hypo_seqs_annotate.fasta Reference_Sequences.faa
```
Diamond will be used for the annotation with the build-in dmnd index.  
The output format mandatory has to be tab separated in order to later script to work.  
The querry sequence ID and the database hit ID in the first and second column respectely.  
```
damiond blastp --db /mnt/f/Cesar_Tesis/MGWAS_Textile/Diamond_uniref90DB/uniref90_diamond.dmnd \
--very-sensitive --outfmt 6 --threads 30 --out Annotated_proteins.tsv \
--query Reference_Sequences.faa
```
* TIP: Perform the annotation over all hypothetical proteins in your database.
* If the Textile Dye database was used the protein annotated file is located at: /mnt/f/Cesar_Tesis/MGWAS_Textile/Metagenoma_Completo/MapeoENS_Bowtie2/Databases/PositiveMG_DB/DB2/Diamond_Annotation/DB2_allproteins.tsv


## Enrich the annotation on the Significant_Gene_Hits.txt file
Finally, we are going to use the Gene_Hits_Annotation_Args.py to combine this annotation on to the Significant_Gene_Hits.txt file.  
We need to execute the script from a fixed location (/mnt/f/Cesar_Tesis/Hits_To_UniRef_KOs/) as there are ten files with .split (approximately 8Gb) extensión that contains unirefIDs and it's annotation.
Also all the KO IDs are also there (KOs.keg). Threfore, the full path for our input and ouput files are required.  
```
cd /mnt/f/Cesar_Tesis/Hits_To_UniRef_KOs/
python Gene_Hits_Annotation_Args.py -i /$working_dir/Significant_Gene_Hits.txt -b /$working_dir/Annotated_proteins.tsv -k KOs.keg -o /$working_dir/Significant_Gene_Hits_Annotated.txt
```
The final output should look somehting like this:  
> TEXDB_S11_29835	2	8.655607726	0.488	0.488	0.2181	UniRef90_A0A3M5PEG9	UniRef90_A0A3M5PEG9 Uncharacterized protein  
> TEXDB_S5_258199	3	9.514278574	0.226	0.226	0.0975	UniRef90_A0A7W8M979	UniRef90_A0A7W8M979 Acyl-CoA dehydrogenase   
> erpA_3	1	11.1739252	0.214	0.214	0.164	K15724  erpA; iron-sulfur cluster insertion protein  	
> ndhH_4	1	8.314258261	0.238	0.238	0.119	K05579  ndhH; NAD(P)H-quinone oxidoreductase subunit H [EC:7.1.1.2]  

The script Add_KOs_UnirefIDs.sh allows to iterate this process over multiple "Significant_Gene_Hits.txt" results that you want to annotated. You only need to create a ouput folder to store the new files
and change the path in the script.
```
mkdir -p Results # This named has to be placed in line 11 of the script before executing
bash Add_KOs_UnirefIDs.sh
```
**THIS CONCLUDE THE EXECUTION OF PYSEER WITH THE LLM FOR BINARY PHENOTYPE**
## De novo constructing a new database reference  
If you want to annotate a new set of metagenomes and used it as refence for Pyseer activate metaprokka enviroment and used the following script once the variables are updated
according to you.  

```
conda activate metaprokka
```
Copy this to a executable file 
```
#!/bin/bash
DATABASE_NAME=your_database_name  # Name of your database 
DATABASE_DIR='path/to/database' # Full path to working directory or where to store the .fasta and .gff of the database
PATTERN_NAME=NEWDB # change to your own
i=1
while read -r line; do
    sample_name=$(basename "$line" .fasta)
    metaprokka --locustag ${PATTERN_NAME}$i --prefix ${PATTERN_NAME}$i --outdir ${sample_name}$i --norrna --cpus 40 "$line"

    if [-d "${sample_name}$i"]; then
        cd "${sample_name}$i"/
        grep "${PATTERN_NAME}$i" "${PATTERN_NAME}$i.gff" >> "$DATABASE_DIR"/"$DATABASE_NAME".gff
        cp "${PATTERN_NAME}$i.ffn" "$DATABASE_DIR"/"${PATTERN_NAME}$i.ffn"
        cd ..
    fi

    i=$((i+1))
done < samples.list #

cd "$DATABASE_DIR/"
cat *.ffn > "$DATABASE_NAME".fasta
rm -f *.ffn
```

# References
Lees, John A., Galardini, M., et al. pyseer: a comprehensive tool for microbial pangenome-wide association studies. 
Bioinformatics 34:4310–4312 (2018). doi:10.1093/bioinformatics/bty539  

Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs. 
bioRxiv 695338 (2019). doi: https://doi.org/10.1101/695338 



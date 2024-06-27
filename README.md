# How-to-perfom-a-MWAS-at-Ayixon-Lab
This repository intens to serve as a blueprint to perform a microbiome-wide association study on the archezoa server.
The main software useded is Pyseer (Lees et al., 2017), how ever others are also used, if not mentioned all are is installed on $PATH.
For further Pyseer information refer to : https://pyseer.readthedocs.io/en/master/tutorial.html 

## Prepare input files

We'll assume that all assemblies are on a folder, hereinafter referred to as "working_folder".
* PRO TIP:
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
## Variant (unitigs) Callin 
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
```
```
## Performing and linear mixed model on Pysser.
```
```
## Variant annotation
```
```
## De novo constructing a database reference
```
```

# References
Lees, John A., Galardini, M., et al. pyseer: a comprehensive tool for microbial pangenome-wide association studies. 
Bioinformatics 34:4310–4312 (2018). doi:10.1093/bioinformatics/bty539 
Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs. 
bioRxiv 695338 (2019). doi: https://doi.org/10.1101/695338 



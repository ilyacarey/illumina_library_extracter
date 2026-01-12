This Nextflow pipeline counts library variants per bin from raw paired Illumina reads. 

**Key steps:**

-fastp for quality control and adapter trimming

-flash2 for merging reads

-cutadapt to orient reads using the inline barcodes and the reverse primer

-cutadapt to demultiplex reads with inline barcodes

-cutadapt to extract a library region using flanking sequences

-a script to count the frequency of each sequence in each bin


The final outputs are .tsv files with sequence frequency counts per bin, as well as a .tsv with all the bins merged. 
A multiqc report is also generated.



**Requirements**:

-java

-nextflow

-conda

To run this script, update the relevant parameters in the config file (or when running in bash), then run in bash by cloning or directly from repo:


git clone https://github.com/ilyacarey/illumina_library_extracter.git

cd REPO

nextflow run main.nf -profile conda




or remotely (in which case you must have inline barcodes fasta file available):

nextflow run ilyacarey/illumina_library_extracter -profile conda 





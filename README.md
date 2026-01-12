This Nextflow pipeline extracts a specific region (e.g., a variant library) from raw paired Illumina reads. 

Key steps:
-fastp for quality control and adapter trimming

-flash2 for merging reads

-cutadapt to orient reads using the inline barcodes and the reverse primer

-cutadapt to demultiplex reads with inline barcodes

-cutadapt to extract a library region using flanking sequences

-a script to count the frequency of each sequence in each bin

The final outputs .tsv files with sequence frequency counts per bin, as well as a .tsv with all the bins merged. 
A multiqc report is also generated.

# illumina_library_extracter

Nextflow pipeline for processing paired-end Illumina reads from sort-seq or similar library screens. Takes raw FASTQ files and outputs per-bin sequence frequency counts.

## Pipeline steps

1. **fastp** — quality control and adapter trimming
2. **FLASH2** — merge paired-end reads by overlap
3. **cutadapt** — orient reads using inline barcodes and a 3′ anchor sequence
4. **cutadapt** — demultiplex into per-bin FASTQ files using inline barcodes
5. **cutadapt** *(optional)* — extract the library region between left and right flanking sequences; reads that are untrimmed, too short, or too long are written to separate files for inspection
6. **Python** — count the frequency of each unique sequence per bin
7. **MultiQC** — aggregate QC report across all steps

## Outputs

- `06_counts/` — per-bin sequence counts for the extracted library region (TSV)
- `07_merged_counts/` — all bins merged into a single TSV
- `06_whole_counts/` — per-bin sequence counts for full demuxed reads (TSV)
- `07_merged_whole_counts/` — all bins merged into a single TSV
- `00_reports/` — MultiQC HTML report

Whole-read and library-region counting can each be toggled independently (see parameters below).

## Requirements

- Java
- Nextflow
- conda

## Usage

**Clone and run locally:**
```bash
git clone https://github.com/ilyacarey/illumina_library_extracter.git
cd illumina_library_extracter
nextflow run main.nf -profile conda
```

**Run directly from GitHub** (barcodes FASTA must be available locally):
```bash
nextflow run ilyacarey/illumina_library_extracter -profile conda
```

Update the relevant parameters in `nextflow.config` before running, or pass them on the command line:
```bash
nextflow run main.nf -profile conda \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --sample my_sample \
  --barcodes_fasta barcodes_anchored.fasta
```

## Key parameters

| Parameter | Description |
|---|---|
| `read1` / `read2` | Input paired-end FASTQ files |
| `sample` | Sample name used for output file naming |
| `barcodes_fasta` | FASTA of inline barcodes (must be 5′-anchored with `^`) |
| `orient_3p_anchor` | 3′ anchor sequence for read orientation (anchor to end with `$`) |
| `left_flank` | Sequence immediately upstream of the library region |
| `right_flank` | Sequence immediately downstream of the library region |
| `lib_min_len` / `lib_max_len` | Expected length range (bp) of the extracted library region |
| `count_library` | Output counts of extracted library region (default: `true`) |
| `count_whole_reads` | Output counts of full demuxed reads (default: `true`) |
| `cleanup` | Delete `work/` directory on successful completion (default: `true`); set to `false` to enable `-resume` |
| `cpus` | Number of CPU threads (default: `8`) |

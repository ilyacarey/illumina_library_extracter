nextflow.enable.dsl=2

/*
 * Required inputs
 */
if( !params.read1 || !params.read2 ) {
  log.error """
  Missing inputs.

  Usage:
    nextflow run main.nf -profile conda \\
      --read1 <R1.fastq.gz> \\
      --read2 <R2.fastq.gz> \\
      --sample <sample_name> \\
      --barcodes_fasta <barcodes_anchored.fasta> \\
      --outdir results
  """
  System.exit(1)
}


/*
 * Single-sample input channel
 */
def r1 = file(params.read1)
def r2 = file(params.read2)
def barcodes_file = file(params.barcodes_fasta)

Channel
  .of( tuple(params.sample, r1, r2) )
  .set { ch_sample }

barcodes_ch = Channel.value(barcodes_file)


/*
 * Processes
 */
process FASTP {
  tag "${sample_id}"
  publishDir "${params.outdir}/01_fastp/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    // Added emit: reads
    tuple val(sample_id), path("${sample_id}_1.fastp.fastq.gz"), path("${sample_id}_2.fastp.fastq.gz"), path("${sample_id}.fastp.html"), path("${sample_id}.fastp.json"), emit: reads
    path("*.json"), emit: log 
    path("*.html")

  script:
  """
  fastp \
    -i ${read1} -I ${read2} \
    -o ${sample_id}_1.fastp.fastq.gz \
    -O ${sample_id}_2.fastp.fastq.gz \
    --detect_adapter_for_pe \
    --cut_front --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --length_required 50 \
    --n_base_limit 5 \
    --thread ${task.cpus} \
    --html ${sample_id}.fastp.html \
    --json ${sample_id}.fastp.json
  """
}

// Renamed process to match workflow call
process FLASH2_MERGE {
  tag "${sample_id}"
  publishDir "${params.outdir}/02_flash2/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    // Input tuple must match the output of FASTP (5 elements) or use just what is needed
    tuple val(sample_id), path(r1), path(r2), path(html), path(json)

  output:
    // Added emit: reads
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), emit: reads
    path "${sample_id}.flash.log", emit: log

  script:
  """
  # Added log redirection (> ... 2>&1) so the log file is actually created
  flash2 \
    -t ${task.cpus} \
    -m ${params.flash_min_overlap} \
    -M ${params.flash_max_overlap} \
    -x ${params.flash_mismatch} \
    -p ${params.flash_phred} \
    -o ${sample_id} \
    ${r1} ${r2} > ${sample_id}.flash.log 2>&1

  gzip -f ${sample_id}.extendedFrags.fastq
  """
}

process ORIENT_READS {
  tag "${sample_id}"
  publishDir "${params.outdir}/03_orient/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(merged)
    path barcodes

  output:
    // Added emit: reads
    tuple val(sample_id), path("${sample_id}.merged.oriented.fastq.gz"), path("${sample_id}.orient.log"), emit: reads
    path "${sample_id}.orient.log", emit: log

  script:
  """
  cutadapt \
    -j ${task.cpus} \
    --revcomp \
    --no-trim \
    -g file:${barcodes} \
    -a ${params.orient_3p_anchor} \
    --discard-untrimmed \
    -o ${sample_id}.merged.oriented.fastq.gz \
    ${merged} \
    > ${sample_id}.orient.log
  """
}

process DEMUX {
  tag "${sample_id}"
  publishDir "${params.outdir}/04_demux/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(oriented), path(orient_log)
    path barcodes

  output:
    // Added emit: reads
    tuple val(sample_id), path("demux/*.fastq.gz"), path("${sample_id}.demux.log"), emit: reads
    path "${sample_id}.demux.log", emit: log

  script:
  """
  mkdir -p demux

  cutadapt \
    -j ${task.cpus} \
    --no-indels \
    --overlap 9 \
    -e 0.1 \
    -g file:${barcodes} \
    -o demux/{name}.fastq.gz \
    --untrimmed-output demux/unknown.fastq.gz \
    ${oriented} \
    > ${sample_id}.demux.log
  """
}

process EXTRACT_TFBS {
  tag "${sample_id}:${bin_fastq.simpleName}"
  publishDir "${params.outdir}/05_extract/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(bin_fastq)

  output:
    tuple val(sample_id),
          path("${bin_fastq.simpleName}.tfbs.fastq.gz"),
          path("${bin_fastq.simpleName}.untrimmed.fastq.gz"),
          path("${bin_fastq.simpleName}.too_short.fastq.gz"),
          path("${bin_fastq.simpleName}.too_long.fastq.gz"),
          path("${bin_fastq.simpleName}.extract.log"),
          emit: reads
    path("${bin_fastq.simpleName}.extract.log"), emit: log

  script:
  """
  cutadapt \
    -j ${task.cpus} \
    --no-indels \
    -e ${params.cutadapt_e} \
    --overlap ${params.extract_overlap} \
    -g "${left}" \
    -a "${right}" \
    -m ${params.tfbs_len} -M ${params.tfbs_len} \
    --discard-untrimmed \
    --untrimmed-output ${bin_fastq.simpleName}.untrimmed.fastq.gz \
    --too-short-output ${bin_fastq.simpleName}.too_short.fastq.gz \
    --too-long-output ${bin_fastq.simpleName}.too_long.fastq.gz \
    -o ${bin_fastq.simpleName}.tfbs.fastq.gz \
    ${bin_fastq} \
    > ${bin_fastq.simpleName}.extract.log
  """
}

}

process COUNT_TFBS {
  tag "${sample_id}:bin${bin_num}"
  publishDir "${params.outdir}/06_counts/${sample_id}", mode: 'copy'
  cpus 1

  input:
    tuple val(sample_id), val(bin_num), path(tfbs_fastq)

  output:
    tuple val(sample_id), path("${tfbs_fastq.simpleName}.counts.tsv")

  script:
  """
  python - << 'PY'
  import gzip
  from collections import Counter

  tfbs_path = "${tfbs_fastq}"
  out_path  = "${tfbs_fastq.simpleName}.counts.tsv"
  bin_num   = "${bin_num}"

  counts = Counter()
  with gzip.open(tfbs_path, "rt") as fh:
    for i, line in enumerate(fh):
      if i % 4 == 1:
        seq = line.strip()
        if seq:
          counts[seq] += 1

  with open(out_path, "w") as out:
    out.write("bin\\tsequence\\tcount\\n")
    for seq, c in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
      out.write(f"{bin_num}\\t{seq}\\t{c}\\n")
  PY
  """
}


process MERGE_COUNTS {
  tag "${sample_id}"
  publishDir "${params.outdir}/07_merged_counts/${sample_id}", mode: 'copy'
  cpus 1

  input:
    tuple val(sample_id), path(count_files)

  output:
    path("all_bins.counts.tsv")

  script:
  // Using a cleaner bash loop to safely handle list of files
  def files_str = count_files.join(' ')
  """
  set -euo pipefail
  
  # Grab header from first file
  first_file=\$(echo "${files_str}" | awk '{print \$1}')
  head -n 1 \$first_file > all_bins.counts.tsv

  # Loop and skip header
  for f in ${files_str}; do
    tail -n +2 "\$f" >> all_bins.counts.tsv
  done
  """
}


process MULTIQC {
    publishDir "${params.outdir}/00_reports", mode: 'copy'

    input:
    path multiqc_files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

/*
 * Workflow wiring
 */
workflow {

  ch_fastp  = FASTP(ch_sample)
  // Use .reads output specifically
  ch_merged = FLASH2_MERGE(ch_fastp.reads)

  // Use .reads output
  ch_oriented = ORIENT_READS(ch_merged.reads, barcodes_ch)
  // Use .reads output
  ch_demux    = DEMUX(ch_oriented.reads, barcodes_ch)

  // Flatten demux bins, skip unknown
  // ch_demux.reads is now the tuple containing the file list
  ch_bins = ch_demux.reads.flatMap { sample_id, files, demux_log ->
    files
      .findAll { it.name != 'unknown.fastq.gz' }
      .collect { f -> tuple(sample_id, f) }
  }

  ch_tfbs = EXTRACT_TFBS(ch_bins)

  // Derive numeric bin from filename
  // Use .reads output
  ch_counts_in = ch_tfbs.reads.map { sample_id, tfbs_fastq, extract_log ->
    def m = (tfbs_fastq.simpleName =~ /(\d+)/)
    def bin = m.find() ? m.group(1) : 'NA'
    tuple(sample_id, bin, tfbs_fastq)
  }

  ch_counts = COUNT_TFBS(ch_counts_in)

  ch_counts
    .groupTuple()
    .map { sample_id, files -> tuple(sample_id, files) }
    | MERGE_COUNTS

  // 1. Gather logs from all steps
  ch_reports = Channel.empty()
      .mix(ch_fastp.log)
      .mix(ch_merged.log)
      .mix(ch_oriented.log)
      .mix(ch_demux.log)
      .mix(ch_tfbs.log)
      .collect() 

  // 2. Generate Report
  MULTIQC(ch_reports)
}

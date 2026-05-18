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

process FLASH2_MERGE {
  tag "${sample_id}"
  publishDir "${params.outdir}/02_flash2/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(r1), path(r2), path(html), path(json)

  output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), emit: reads
    path "${sample_id}.flash.log", emit: log

  script:
  """
  flash2 \
    -z \
    -t ${task.cpus} \
    -m ${params.flash_min_overlap} \
    -M ${params.flash_max_overlap} \
    -x ${params.flash_mismatch} \
    -p ${params.flash_phred} \
    -o ${sample_id} \
    ${r1} ${r2} > ${sample_id}.flash.log 2>&1
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
    tuple val(sample_id), path("${sample_id}.merged.oriented.fastq.gz"), path("${sample_id}.orient.log"), emit: reads
    path "${sample_id}.orient.log", emit: log

  script:
  """
  cutadapt \
    -j ${task.cpus} \
    --revcomp \
    --no-trim \
    -g file:${barcodes} \
    -e ${params.cutadapt_e_orient} \
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
      tuple val(sample_id), path("demux/*.fastq.gz"), path("${sample_id}_demux.cutadapt.log"), emit: reads
      path "${sample_id}_demux.cutadapt.log", emit: log

  script:
  """
  mkdir -p demux

  cutadapt \
    -j ${task.cpus} \
    --overlap 9 \
    -e ${params.cutadapt_e_demux} \
    -g file:${barcodes} \
    -o demux/{name}.fastq.gz \
    --untrimmed-output demux/unknown.fastq.gz \
    ${oriented} \
    > ${sample_id}_demux.cutadapt.log
  """
}

process EXTRACT_LIBRARY {
  tag "${sample_id}:${bin_fastq.simpleName}"
  publishDir "${params.outdir}/05_extract/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(bin_fastq)

  output:
    tuple val(sample_id), path("${bin_fastq.simpleName}.library.fastq.gz"), path("${bin_fastq.simpleName}.extract.log"), emit: reads
    path "*.extract.log", emit: log
    path "${bin_fastq.simpleName}.untrimmed.fastq.gz",  optional: true
    path "${bin_fastq.simpleName}.toolong.fastq.gz",    optional: true
    path "${bin_fastq.simpleName}.tooshort.fastq.gz",   optional: true

  script:
  """
  cutadapt \
    -j ${task.cpus} \
    --no-indels \
    -e ${params.cutadapt_e} \
    --overlap ${params.extract_overlap} \
    -g "${params.left_flank}...${params.right_flank}" \
    --untrimmed-output ${bin_fastq.simpleName}.untrimmed.fastq.gz \
    -m ${params.lib_min_len} -M ${params.lib_max_len} \
    --too-long-output  ${bin_fastq.simpleName}.toolong.fastq.gz \
    --too-short-output ${bin_fastq.simpleName}.tooshort.fastq.gz \
    -o ${bin_fastq.simpleName}.library.fastq.gz \
    ${bin_fastq} \
    > ${bin_fastq.simpleName}.extract.log
  """
}

process COUNT_LIBRARY {
  tag "${sample_id}:bin${bin_num}"
  publishDir "${params.outdir}/06_counts/${sample_id}", mode: 'copy'
  cpus 1

  input:
    tuple val(sample_id), val(bin_num), path(lib_fastq)

  output:
    tuple val(sample_id), path("${lib_fastq.simpleName}.counts.tsv")

  script:
  """
  python - << 'PY'
  import gzip
  from collections import Counter

  lib_path = "${lib_fastq}"
  out_path  = "${lib_fastq.simpleName}.counts.tsv"
  bin_num   = "${bin_num}"

  counts = Counter()
  with gzip.open(lib_path, "rt") as fh:
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


process COUNT_WHOLE_READS {
  tag "${sample_id}:bin${bin_num}"
  publishDir "${params.outdir}/06_whole_counts/${sample_id}", mode: 'copy'
  cpus 1

  input:
    tuple val(sample_id), val(bin_num), path(fastq)

  output:
    tuple val(sample_id), path("${fastq.simpleName}.counts.tsv")

  script:
  """
  python - << 'PY'
  import gzip
  from collections import Counter

  fastq_path = "${fastq}"
  out_path   = "${fastq.simpleName}.counts.tsv"
  bin_num    = "${bin_num}"

  counts = Counter()
  with gzip.open(fastq_path, "rt") as fh:
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


process MERGE_WHOLE_COUNTS {
  tag "${sample_id}"
  publishDir "${params.outdir}/07_merged_whole_counts/${sample_id}", mode: 'copy'
  cpus 1

  input:
    tuple val(sample_id), path(count_files)

  output:
    path("all_bins.counts.tsv")

  script:
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
  ch_merged = FLASH2_MERGE(ch_fastp.reads)

  ch_oriented = ORIENT_READS(ch_merged.reads, barcodes_ch)
  ch_demux    = DEMUX(ch_oriented.reads, barcodes_ch)

  // Flatten demux bins, skip unknown
  ch_bins = ch_demux.reads.flatMap { sample_id, files, demux_log ->
    files
      .findAll { it.name != 'unknown.fastq.gz' }
      .collect { f -> tuple(sample_id, f) }
  }

  // Count whole (demuxed) reads per bin
  if (params.count_whole_reads) {
    ch_bins_numbered = ch_bins.map { sample_id, bin_fastq ->
      def m = (bin_fastq.simpleName =~ /(\d+)/)
      def bin = m.find() ? m.group(1) : 'NA'
      tuple(sample_id, bin, bin_fastq)
    }
    COUNT_WHOLE_READS(ch_bins_numbered)
      .groupTuple()
      .map { sample_id, files -> tuple(sample_id, files) }
      | MERGE_WHOLE_COUNTS
  }

  // Gather base logs
  ch_reports = Channel.empty()
      .mix(ch_fastp.log)
      .mix(ch_merged.log)
      .mix(ch_oriented.log)
      .mix(ch_demux.log)

  // Extract library region and count per bin
  if (params.count_library) {
    ch_lib = EXTRACT_LIBRARY(ch_bins)

    ch_lib.reads
      .map { sample_id, lib_fastq, extract_log ->
        def m = (lib_fastq.simpleName =~ /(\d+)/)
        def bin = m.find() ? m.group(1) : 'NA'
        tuple(sample_id, bin, lib_fastq)
      }
      | COUNT_LIBRARY
      | groupTuple
      | map { sample_id, files -> tuple(sample_id, files) }
      | MERGE_COUNTS

    ch_reports = ch_reports.mix(ch_lib.log)
  }

  // Generate Report
  MULTIQC(ch_reports.collect())
}

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

  Example:
    nextflow run main.nf -profile conda \\
      --read1 NG-A3433_CarR_de_novo_libLAJ7476_1.fastq.gz \\
      --read2 NG-A3433_CarR_de_novo_libLAJ7476_2.fastq.gz \\
      --sample NG-A3433 \\
      --barcodes_fasta barcodes_anchored.fasta \\
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
    tuple val(sample_id),
      path("${sample_id}_1.fastp.fastq.gz"),
      path("${sample_id}_2.fastp.fastq.gz"),
      path("${sample_id}.fastp.html"),
      path("${sample_id}.fastp.json")

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

process FLASH_MERGE {
  tag "${sample_id}"
  publishDir "${params.outdir}/02_flash/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(r1), path(r2), path(html), path(json)

  output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz")

  script:
  """
  mkdir -p flash_out

  flash \
    -t ${task.cpus} \
    -m ${params.flash_min_overlap} \
    -M ${params.flash_max_overlap} \
    -x ${params.flash_mismatch} \
    -p ${params.flash_phred} \
    -o ${sample_id} \
    -d flash_out \
    ${r1} ${r2}

  gzip -f flash_out/${sample_id}.extendedFrags.fastq
  mv flash_out/${sample_id}.extendedFrags.fastq.gz ${sample_id}.extendedFrags.fastq.gz
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
    tuple val(sample_id), path("${sample_id}.merged.oriented.fastq.gz"), path("${sample_id}.orient.log")

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
    tuple val(sample_id), path("demux/*.fastq.gz"), path("${sample_id}.demux.log")

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
    tuple val(sample_id), path("${bin_fastq.simpleName}.tfbs.fastq.gz"), path("${bin_fastq.simpleName}.extract.log")

  script:
  """
  cutadapt \
    -j ${task.cpus} \
    --no-indels \
    -e ${params.cutadapt_e} \
    --overlap ${params.extract_overlap} \
    -g "${params.left_flank}...${params.right_flank}" \
    --discard-untrimmed \
    -m ${params.tfbs_len} -M ${params.tfbs_len} \
    -o ${bin_fastq.simpleName}.tfbs.fastq.gz \
    ${bin_fastq} \
    > ${bin_fastq.simpleName}.extract.log
  """
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
      # FASTQ sequence line is every 4 lines: 2nd line in each record (i % 4 == 1)
      if i % 4 == 1:
        seq = line.strip()
        if seq:
          counts[seq] += 1

  with open(out_path, "w") as out:
    out.write("bin\\tsequence\\tcount\\n")
    # sort by count descending, then sequence (stable, nice for diffs)
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
  """
  set -euo pipefail

  head -n 1 ${count_files[0]} > all_bins.counts.tsv
  for f in ${count_files}; do
    tail -n +2 "\$f" >> all_bins.counts.tsv
  done
  """
}


/*
 * Workflow wiring
 */
workflow {

  ch_fastp  = FASTP(ch_sample)
  ch_merged = FLASH_MERGE(ch_fastp)

  ch_oriented = ORIENT_READS(ch_merged, barcodes_ch)
  ch_demux    = DEMUX(ch_oriented, barcodes_ch)

  // Flatten demux bins, skip unknown
  ch_bins = ch_demux.flatMap { sample_id, files, demux_log ->
    files
      .findAll { it.name != 'unknown.fastq.gz' }
      .collect { f -> tuple(sample_id, f) }
  }

  ch_tfbs = EXTRACT_TFBS(ch_bins)

  // Derive numeric bin from filename in Groovy (avoids bash parsing problems)
  ch_counts_in = ch_tfbs.map { sample_id, tfbs_fastq, extract_log ->
    def m = (tfbs_fastq.simpleName =~ /(\d+)/)
    def bin = m.find() ? m.group(1) : 'NA'
    tuple(sample_id, bin, tfbs_fastq)
  }

  ch_counts = COUNT_TFBS(ch_counts_in)

  ch_counts
    .groupTuple()
    .map { sample_id, files -> tuple(sample_id, files) }
    | MERGE_COUNTS
}

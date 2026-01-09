nextflow.enable.dsl=2

/*
 * Default parameters
 */
params.reads_dir        = params.reads_dir ?: 'data'
params.reads_pattern    = params.reads_pattern ?: '*_{1,2}.fastq.gz'
params.outdir           = params.outdir ?: 'results'

params.cpus             = params.cpus ?: 8

params.barcodes_fasta   = params.barcodes_fasta ?: 'barcodes_anchored.fasta'

params.orient_3p_anchor = params.orient_3p_anchor ?: 'GATCTATGAGCAAAGGAGAAGAAC$'

params.left_flank       = params.left_flank ?: 'TTTACGGCTAGCTCAGTCCTAGGTACAATGCTAGCGAATTC'
params.right_flank      = params.right_flank ?: 'GGATCCAGC'

params.cutadapt_e       = params.cutadapt_e ?: 0.05
params.extract_overlap  = params.extract_overlap ?: 9
params.tfbs_len         = params.tfbs_len ?: 20

// Your FLASH settings (keep as params so you can tune later)
params.flash_min_overlap = params.flash_min_overlap ?: 20
params.flash_max_overlap = params.flash_max_overlap ?: 150
params.flash_mismatch    = params.flash_mismatch ?: 0.1
params.flash_phred       = params.flash_phred ?: 33


/*
 * Input detection: auto-pair *_1.fastq.gz with *_2.fastq.gz
 */
Channel
  .fromFilePairs("${params.reads_dir}/${params.reads_pattern}", flat: true)
  .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
  .set { ch_samples }

// Stage barcode fasta into work dirs reliably
barcodes_ch = Channel.value( file(params.barcodes_fasta) )

/*
 * Processes
 */
process FASTP {
  tag "${sample_id}"
  publishDir "${params.outdir}/fastp", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id),
          path("${sample_id}.R1.fastp.fq.gz"),
          path("${sample_id}.R2.fastp.fq.gz"),
          path("${sample_id}.fastp.html"),
          path("${sample_id}.fastp.json")

  script:
  """
  fastp \
    -i ${r1} -I ${r2} \
    -o ${sample_id}.R1.fastp.fq.gz \
    -O ${sample_id}.R2.fastp.fq.gz \
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
  publishDir "${params.outdir}/flash", mode: 'copy'
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
  publishDir "${params.outdir}/orient", mode: 'copy'
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
  publishDir "${params.outdir}/demux/${sample_id}", mode: 'copy'
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
  tag "${sample_id}:${bin_file.simpleName}"
  publishDir "${params.outdir}/extract/${sample_id}", mode: 'copy'
  cpus params.cpus

  input:
    tuple val(sample_id), path(bin_file)

  output:
    tuple val(sample_id),
          path("${bin_file.simpleName}.tfbs.fastq.gz"),
          path("${bin_file.simpleName}.extract.log")

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
    -o ${bin_file.simpleName}.tfbs.fastq.gz \
    ${bin_file} \
    > ${bin_file.simpleName}.extract.log
  """
}

process COUNT_TFBS {
  tag "${sample_id}:${tfbs_fastq.simpleName}"
  publishDir "${params.outdir}/counts/${sample_id}", mode: 'copy'
  cpus 1

  input:
    tuple val(sample_id), path(tfbs_fastq)

  output:
    tuple val(sample_id), path("${tfbs_fastq.simpleName}.counts.tsv")

  script:
  """
  set -euo pipefail

  bin_name="${tfbs_fastq.simpleName}"
  bin_num="NA"
  if [[ "\$bin_name" =~ ([0-9]+) ]]; then
    bin_num="\${BASH_REMATCH[1]}"
  fi

  {
    echo -e "bin\\tsequence\\tcount"
    gzip -dc ${tfbs_fastq} \
      | awk 'NR%4==2' \
      | sort \
      | uniq -c \
      | awk -v b="\$bin_num" '{print b "\\t" \$2 "\\t" \$1}'
  } > ${tfbs_fastq.simpleName}.counts.tsv
  """
}

process MERGE_COUNTS {
  tag "${sample_id}"
  publishDir "${params.outdir}/counts/${sample_id}", mode: 'copy'
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

  ch_fastp   = FASTP(ch_samples)
  ch_merged  = FLASH_MERGE(ch_fastp)

  // pass barcode file into processes that need it
  ch_oriented = ORIENT_READS(ch_merged, barcodes_ch)
  ch_demux    = DEMUX(ch_oriented, barcodes_ch)

  // Flatten: (sample_id, [many demux fastqs]) -> (sample_id, each_bin_fastq)
  ch_bins = ch_demux
    .flatMap { sample_id, files, demux_log ->
      files
        .findAll { it.name != 'unknown.fastq.gz' }
        .collect { f -> tuple(sample_id, f) }
    }

  ch_tfbs = EXTRACT_TFBS(ch_bins)

  ch_counts = COUNT_TFBS(ch_tfbs.map { sample_id, tfbs_fastq, log -> tuple(sample_id, tfbs_fastq) })

  ch_counts
    .groupTuple()
    .map { sample_id, files -> tuple(sample_id, files) }
    | MERGE_COUNTS
}

nextflow.enable.dsl=2

// Required inputs
if( !params.read1 || !params.read2 ) {
  log.error "Missing inputs. Usage:\n" +
            "  nextflow run main.nf -profile conda --read1 <R1.fastq.gz> --read2 <R2.fastq.gz> --outdir results\n"
  System.exit(1)
}

params.outdir = params.outdir ?: 'results'
params.cpus   = (params.cpus ?: 8) as int
params.sample = params.sample ?: 'sample'

// Convert to Nextflow file objects
def r1 = file(params.read1)
def r2 = file(params.read2)

Channel
  .of( tuple(params.sample, r1, r2) )
  .set { ch_sample }

process FASTP {
  tag "${sample_id}"
  publishDir "${params.outdir}/fastp/${sample_id}", mode: 'copy'
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

workflow {
  FASTP(ch_sample)
}

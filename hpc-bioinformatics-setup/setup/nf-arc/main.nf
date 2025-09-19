nextflow.enable.dsl=2

// ---------------- Params ----------------
params.sample_id  = params.sample_id ?: 'sample1'
params.libraries  = params.libraries ?: "${System.getenv('LAB_SCRATCH')}/projects/multiome/${params.sample_id}/libraries.csv"
params.reference  = params.reference ?: "${System.getenv('LAB_ARCHIVE')}/ref/refdata-cellranger-arc-GRCh38-2020-A"
params.outdir     = params.outdir    ?: "${System.getenv('LAB_SCRATCH')}/results/arc/${params.sample_id}"
params.extra_args = params.extra_args ?: ''   // place for --expect-cells etc.

// ---------------- Workflow --------------
workflow {
  RUN_ARC_COUNT(params.sample_id, file(params.libraries), params.reference)
}

// ---------------- Process ----------------
process RUN_ARC_COUNT {
  tag "${sample_id}"
  cpus 16
  memory '64 GB'
  time '24h'
  // Cell Ranger does its own threading; map NF resources to CR args:
  // publish the essential outputs
  publishDir params.outdir, mode: 'copy', createDir: true, overwrite: true

  input:
  val sample_id
  path libs_csv
  val ref_path

  output:
  // copy the whole outs/ to result dir; but also explicitly expose key files
  path "${sample_id}/outs", emit: outs_dir
  path "${sample_id}/outs/atac_fragments.tsv.gz", optional: true, emit: atac_frags
  path "${sample_id}/outs/atac_fragments.tsv.gz.tbi", optional: true
  path "${sample_id}/outs/peaks.bed", optional: true, emit: atac_peaks
  path "${sample_id}/outs/filtered_feature_bc_matrix.h5", optional: true, emit: gex_h5
  path "${sample_id}/outs/filtered_feature_bc_matrix", optional: true
  path "${sample_id}/outs/summary.csv", optional: true

  script:
  // translate NF resource specs to what cellranger-arc expects
  def cores = task.cpus
  def memGB = task.memory.toGiga()
  """
  set -euo pipefail

  export TMPDIR="${System.getenv('TMPDIR') ?: System.getenv('LAB_SCRATCH') + '/tmp'}"

  # ensure we see your tools dir from ~/.bashrc; fallback if not sourced
  export PATH="$HOME/tools/cellranger-arc-2.0.2:\$PATH"

  cellranger-arc count \
    --id=${sample_id} \
    --libraries=${libs_csv} \
    --reference=${ref_path} \
    --localcores=${cores} \
    --localmem=${memGB} \
    ${params.extra_args}

  # publish key files to the top of outdir for convenience
  mkdir -p publish_tmp
  cp -f ${sample_id}/outs/summary.csv publish_tmp/ 2>/dev/null || true
  cp -f ${sample_id}/outs/atac_fragments.tsv.gz* publish_tmp/ 2>/dev/null || true
  cp -f ${sample_id}/outs/peaks.bed publish_tmp/ 2>/dev/null || true
  cp -f ${sample_id}/outs/filtered_feature_bc_matrix.h5 publish_tmp/ 2>/dev/null || true

  """
}

workflow.onComplete {
  println "Done. Results in: ${params.outdir}"
  println "  - outs/ directory with full Cell Ranger ARC outputs"
  println "  - convenience copies: fragments.tsv.gz(.tbi), peaks.bed, filtered_feature_bc_matrix.h5, summary.csv"
}

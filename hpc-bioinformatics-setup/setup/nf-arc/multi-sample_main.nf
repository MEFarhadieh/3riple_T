nextflow.enable.dsl=2

// ---------------- Params ----------------
params.sample_id   = params.sample_id   ?: 'sample1'    // used only for single-sample mode
params.libraries   = params.libraries   ?: "${System.getenv('LAB_SCRATCH')}/projects/multiome/${params.sample_id}/libraries.csv"
params.reference   = params.reference   ?: "${System.getenv('LAB_ARCHIVE')}/ref/refdata-cellranger-arc-GRCh38-2020-A"
params.outdir_root = params.outdir_root ?: "${System.getenv('LAB_SCRATCH') ?: baseDir}/results/arc"
// If provided, enables multi-sample mode. CSV must have header: sample_id,libraries
params.samplesheet = params.samplesheet ?: null
params.extra_args  = params.extra_args  ?: ''   // e.g., '--expect-cells=8000'

// ---------------- Channels ----------------
Channel
  .value(params.samplesheet)
  .map { ss ->
      if (ss) {
        // MULTI-SAMPLE MODE from samplesheet
        Channel
          .fromPath(ss)
          .splitCsv(header: true)
          .map { row ->
            def sid   = (row.sample_id ?: row.sample ?: 'NA').toString()
            def libs  = file(row.libraries)
            def refp  = (row.reference ?: params.reference).toString()
            def outd  = "${params.outdir_root}/${sid}"
            tuple(sid, libs, refp, outd)
          }
      } else {
        // SINGLE-SAMPLE MODE (backward compatible)
        def sid  = params.sample_id
        def libs = file(params.libraries)
        def refp = params.reference
        def outd = "${params.outdir_root}/${sid}"
        Channel.of( tuple(sid, libs, refp, outd) )
      }
  }
  .flatten()
  .set { SAMPLES }

// ---------------- Workflow --------------
workflow {
  SAMPLES | RUN_ARC_COUNT
}

// ---------------- Process ----------------
process RUN_ARC_COUNT {
  tag "${sample_id}"
  cpus 16
  memory '64 GB'
  time '24h'

  // publish each sample to its own folder under outdir_root
  publishDir "${outdir}", mode: 'copy', createDir: true, overwrite: true

  input:
  tuple val(sample_id), path(libs_csv), val(ref_path), val(outdir)

  output:
  path "${sample_id}/outs", emit: outs_dir
  path "${sample_id}/outs/atac_fragments.tsv.gz", optional: true, emit: atac_frags
  path "${sample_id}/outs/atac_fragments.tsv.gz.tbi", optional: true
  path "${sample_id}/outs/peaks.bed", optional: true, emit: atac_peaks
  path "${sample_id}/outs/filtered_feature_bc_matrix.h5", optional: true, emit: gex_h5
  path "${sample_id}/outs/filtered_feature_bc_matrix", optional: true
  path "${sample_id}/outs/summary.csv", optional: true

  script:
  def cores = task.cpus
  def memGB = task.memory.toGiga()
  """
  set -euo pipefail

  export TMPDIR="${System.getenv('TMPDIR') ?: System.getenv('LAB_SCRATCH') + '/tmp'}"

  # Prefer module if your HPC uses Environment Modules; fallback to PATH
  if command -v module >/dev/null 2>&1; then
    module purge || true
    module load cellranger-arc/2.0.2 || true
  fi
  export PATH="$HOME/tools/cellranger-arc-2.0.2:\$PATH"

  cellranger-arc count \
    --id=${sample_id} \
    --libraries=${libs_csv} \
    --reference=${ref_path} \
    --localcores=${cores} \
    --localmem=${memGB} \
    ${params.extra_args}

  # convenience copies at the top of the sample's outdir
  mkdir -p publish_tmp
  cp -f ${sample_id}/outs/summary.csv publish_tmp/ 2>/dev/null || true
  cp -f ${sample_id}/outs/atac_fragments.tsv.gz* publish_tmp/ 2>/dev/null || true
  cp -f ${sample_id}/outs/peaks.bed publish_tmp/ 2>/dev/null || true
  cp -f ${sample_id}/outs/filtered_feature_bc_matrix.h5 publish_tmp/ 2>/dev/null || true
  """
}

workflow.onComplete {
  println "Done."
  println "Per-sample results under: ${params.outdir_root}/<sample_id>"
}

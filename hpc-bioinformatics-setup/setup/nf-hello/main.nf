nextflow.enable.dsl=2

params.outdir = params.outdir ?: "${System.getenv('LAB_SCRATCH') ?: baseDir}/results/hello-nf"
def OUTDIR = params.outdir

// ---------- workflow wiring ----------
workflow {
  names = Channel.of('world')     // make the input channel
  SAY_HELLO(names)                
}

// ---------- process ----------
process SAY_HELLO {
  cpus 1
  memory '512 MB'
  time '5m'
  publishDir OUTDIR, mode: 'copy', createDir: true, overwrite: true

  input:
  val name

  output:
  path "hello_${name}.txt"

  script:
  """
  echo "Hello, ${name}!" > hello_${name}.txt
  """
}

workflow.onComplete {
  println "Done. See: ${OUTDIR}"
}

nextflow.enable.dsl=2

def OUTDIR = params.outdir ?: "${System.getenv('LAB_SCRATCH') ?: "$baseDir"}/results/hello-nf"

// tiny input channel
Channel.of('world').set { names }

process SAY_HELLO {
    cpus 1
    memory '512 MB'
    time '5m'
    publishDir OUTDIR, mode: 'copy', overwrite: true

    input:
    val name from names

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

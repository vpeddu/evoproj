process liftOver { 
//conda "${baseDir}/env/env.yml"
//publishDir "${params.OUTPUT}/fastp_PE/${base}", mode: 'symlink', overwrite: true
container "quay.io/biocontainers/ucsc-liftover:377--ha8a8165_4"
beforeScript 'chmod o+rw .'
input: 
    tuple val(species), file(chain_from_human), file(chain_to_human), val(order)
    file human_bed

output: 
    tuple val(species), file("*.check.lifted.bed"), file("*.check.unlifted.bed")

script:
"""
#!/bin/bash

ls -lah

/usr/local/bin/liftOver ${human_bed} ${chain_from_human} ${species}.first.lifted.bed ${species}.first.unlifted.bed

/usr/local/bin/liftOver ${species}.first.lifted.bed ${chain_to_human} ${species}.check.lifted.bed ${species}.check.unlifted.bed

"""
}
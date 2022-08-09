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
    file ("*.check.lifted.bed")

script:
"""
#!/bin/bash

ls -lah

/usr/local/bin/liftOver ${human_bed} ${chain_from_human} ${species}.first.lifted.bed ${species}.first.unlifted.bed

/usr/local/bin/liftOver ${species}.first.lifted.bed ${chain_to_human} ${species}.check.lifted.bed ${species}.check.unlifted.bed

"""
}

process tebag_intersect { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/TEbAG_DB/", mode: 'copy', overwrite: true
container "vpeddu/tebag:latest"
beforeScript 'chmod o+rw .'
input: 
    file lifted_beds
    file human_bed
    file intersect_elements_script

output: 
    tuple file("TEbag_DB.pkl"), file("intersection_table.csv")

script:
"""
#!/bin/bash

ls -lah

python3 ${intersect_elements_script} ${human_bed}

"""
}
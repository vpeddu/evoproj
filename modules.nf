process LiftOver { 
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

process Evo_intersect { 
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

process Evo_match { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/TEbAG_match/", mode: 'copy', overwrite: true
container "vpeddu/tebag:latest"
beforeScript 'chmod o+rw .'
input: 
    file database_pkl
    file quant_file
    file match_script

output: 
    file 'merged_species_TE_df.TEBAG.csv'

script:
"""
#!/bin/bash

ls -lah

python3 ${match_script} ${database_pkl} ${quant_file}
"""
}

process Upset_plot { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/upset_plot/", mode: 'copy', overwrite: true
container "vpeddu/tebag:latest"
beforeScript 'chmod o+rw .'
input: 
    file upset_plot_script
    file TEBAG_MERGED

output: 
    file 'upset_plot.pdf'

"""
#!/bin/bash

ls -lah

Rscript --vanilla ${upset_plot_script} ${TEBAG_MERGED}
"""
}

process Hal2chain_to_human { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/upset_plot/", mode: 'copy', overwrite: true
container "quay.io/comparative-genomics-toolkit/cactus:v2.6.7"
beforeScript 'chmod o+rw .'
cpus 16 
memory '32 GB'
input: 
    val name
    file halfile
    val human_name

output: 
    tuple val(name), file("${name}.human.output/${name}.human.chain.gz")
"""
#!/bin/bash
ls -lah

    cactus-hal2chains --refGenome ${name} \
        --targetGenomes ${human_name} \
        --defaultCores ${task.cpus} \
        --defaultMemory 24G \
        --maxMemory 32G \
        --latest chain_tmp ${halfile} ${name}.human.output
    mv ${name}.human.output/${human_name}.chain.gz ${name}.human.output/${name}.human.chain.gz
"""
}

process Hal2chain_from_human { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/upset_plot/", mode: 'copy', overwrite: true
container "quay.io/comparative-genomics-toolkit/cactus:v2.6.7"
beforeScript 'chmod o+rw .'
cpus 16 
memory '32 GB'
input: 
    val name
    file halfile
    val human_name

output: 
    tuple val(name), file("human.${name}.output/human.${name}.chain.gz")

"""
#!/bin/bash
ls -lah

    cactus-hal2chains --refGenome ${human_name} \
        --targetGenomes ${name} \
        --defaultCores ${task.cpus} \
        --defaultMemory 24G \
        --maxMemory 32G \
        --latest chain_tmp ${halfile} human.${name}.output
    mv human.${name}.output/${name}.chain.gz human.${name}.output/human.${name}.chain.gz
"""
}

process LiftOver_hal { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/LiftOver/${species}", mode: 'copy', overwrite: true
container "quay.io/biocontainers/ucsc-liftover:377--ha8a8165_4"
beforeScript 'chmod o+rw .'
input: 
    tuple val(species), file(chain_from_human), file(chain_to_human)
    file human_bed

output: 
    tuple val(species), file("*.check.lifted.bed"), file("*.check.unlifted.bed")
    file ("*.check*.bed")
    file ("*.first*.bed")

script:
"""
#!/bin/bash

ls -lah

echo lifting ${species}

/usr/local/bin/liftOver ${human_bed} \
    *.human.chain.gz \
    ${species}.first.lifted.bed \
    ${species}.first.unlifted.bed

/usr/local/bin/liftOver ${species}.first.lifted.bed \
    human.*.chain.gz \
    ${species}.check.lifted.bed \
    ${species}.check.unlifted.bed

"""
}

process Splitbed { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/permutation_testing/bed_split/", mode: 'copy', overwrite: true
container "vpeddu/evoproj:v1"
beforeScript 'chmod o+rw .'
input: 
    file annotated_bed

output: 
    path("*.bed", emit: split_bed)
    //path("under_threshold", emit: under_threshold)
shell:
'''
    #!/bin/bash
    ls -lah

    # Process bed file with a simpler approach
    mkdir -p under_threshold
    
 awk '
BEGIN { FS="\\t" }  # Set field separator to tab
$4 ~ /_/ {
    # Extract repeat (first part before underscore)
    repeat = substr($4, 1, index($4, "_")-1)
    
    # Extract species (part after last underscore)
    species_start = 0
    for(i=1; i<=length($4); i++) {
        if(substr($4, i, 1) == "_") 
            species_start = i
    }
    species = substr($4, species_start+1)
    
    # Remove any extensions from species
    if(index(species, ".") > 0)
        species = substr(species, 1, index(species, ".")-1)
    
    # Write to appropriate file
    print $0 > repeat "_" species ".bed"
    
    # Keep track of counts
    counts[repeat "_" species]++
}' !{annotated_bed}

rm *\\)n*

# Move files under threshold (using a single command)
for file in *.bed; do
    if [ -f "$file" ]; then
        # Use stat to get the file size quickly or wc -l
        count=$(awk 'END {print NR}' "$file")
        if [ "$count" -lt 50 ]; then
            mv "$file" "under_threshold/${file%.bed}_under_threshold.bed"
        fi
    fi
done
'''
}
process Runpermutation { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.OUTPUT}/permutation_testing/results/", mode: 'copy', overwrite: true
container "vpeddu/evoproj:v1"
beforeScript 'chmod o+rw .'

cpus 7
memory '16 GB'

input: 
    tuple path(annotated_bed), val(tf), path(bigbeds)
    file fai

output:
    path("*_permutation_results.csv")

script:
"""
#!/bin/bash

    ls -lah
    base=\$(basename "${annotated_bed}")
    IFS='_' read -r -a parts <<< "\$base"
    repeat="\${parts[0]}"
    species="\${parts[-1]%.bed}"

python3 ${baseDir}/bin/permutation_test.py \
    "${annotated_bed}" \
    . \
    ${fai} \
    "\${repeat}" \
    "\${species}" \
    "${tf}" \
    10000
"""
}

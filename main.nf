#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
TEbAG: Transposable Elements by Age Group 
Usage:
quick run command: 
  nextflow run vpeddu/TEBAG \
		 -with-docker 'ubuntu:18.04' \
		 -with-report \
		 -latest \
		 -resume
    """.stripIndent()
}

// show help message
params.help = false
// The params scope allows you to define parameters that will be accessible in the pipeline script
if (params.help){
    helpMessage()
    exit 0
}

include { liftOver } from './modules.nf'

    workflow{
        if ( params.generate_db ){
        Generate_ch = Channel
            .fromPath(params.species_paths)
            // can't get parser to work with headers
            // workaround for now 
            .splitCsv(header: false, skip:1)
            .map { row -> [row[0], file(row[1]), file(row[2]), row[3]] }

        Human_bed = file(params.human_bed)

        liftOver( 
            Generate_ch,
            Human_bed
        )

        }
    }
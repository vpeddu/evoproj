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

// TODO: make modules file specific for tebag generate
include { liftOver } from './modules.nf'
include { tebag_intersect } from './modules.nf'
include { tebag_match} from './modules.nf'

params.generate_db = False

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
        tebag_intersect( 
            liftOver.out[1].collect(),
            Human_bed,
            file("${baseDir}/bin/intersect_elements.py")


        )
        }
        else{ 

            tebag_match(
            file(params.tebag_db)
            file(params.quantification_file)
            file("${baseDir}/bin/tebag_match.py")
            )

        }
    }
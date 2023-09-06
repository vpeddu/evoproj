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
include { LiftOver } from './modules.nf'
include { Tebag_intersect } from './modules.nf'
include { Tebag_match} from './modules.nf'
include { Upset_plot } from './modules.nf'
include { Hal2chain_to_human } from './modules.nf'
include { Hal2chain_from_human } from './modules.nf'
include { LiftOver_hal } from './modules.nf'

params.generate_db = false

    workflow{
        if ( params.generate_db ){
            if (params.CACTUS) { 

            Species_name_Ch = Channel
                .fromPath( params.NAMES_CSV )
                .splitText()
            
            Human_bed = file(params.human_bed)

            Hal2chain_to_human(
                Species_name_Ch,
                file(params.HAL_FILE)
            )
            Hal2chain_from_human(
                Species_name_Ch,
                file(params.HAL_FILE)
            )

            LiftOver_hal(
                Hal2chain_to_human.out.mix(Hal2chain_from_human.out).groupTuple(size:2),
                Human_bed
            )
            Tebag_intersect( 
                LiftOver_hal.out[1].collect(),
                Human_bed,
                file("${baseDir}/bin/intersect_elements.py")
            )

            }
        else { 
            Generate_ch = Channel
                .fromPath(params.species_paths)
                // can't get parser to work with headers
                // workaround for now 
                .splitCsv(header: false, skip:1)
                .map { row -> [row[0], file(row[1]), file(row[2]), row[3]] }

            Human_bed = file(params.human_bed)

            LiftOver( 
                Generate_ch,
                Human_bed
            )
            Tebag_intersect( 
                LiftOver.out[1].collect(),
                Human_bed,
                file("${baseDir}/bin/intersect_elements.py")
            )
            }
        }
        else{ 
            Tebag_match(
            file(params.tebag_db),
            file(params.quantification_file),
            file("${baseDir}/bin/tebag_match.py")
            )
            Upset_plot(
            file("${baseDir}/bin/create_upset_plot.R"),
            Tebag_match.out
            )
        }
    }
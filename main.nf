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
include { Evo_intersect } from './modules.nf'
include { Evo_match} from './modules.nf'
include { Upset_plot } from './modules.nf'
include { Hal2chain_to_human } from './modules.nf'
include { Hal2chain_from_human } from './modules.nf'
include { LiftOver_hal } from './modules.nf'
include {Splitbed} from './modules.nf'
include {Runpermutation} from './modules.nf'

params.generate_db = false
params.CHAINS = false
params.tebag_db = ''

    workflow{
        if ( params.generate_db ){
            if (params.CACTUS) { 

            Species_name_Ch = Channel
                .fromPath( params.NAMES_CSV )
                .splitText().map { it.trim() }
            
            Human_bed = file(params.human_bed)

            Hal2chain_to_human(
                Species_name_Ch,
                file(params.HAL_FILE),
                params.HUMAN_NAME
            )
            Hal2chain_from_human(
                Species_name_Ch,
                file(params.HAL_FILE),
                params.HUMAN_NAME
            )

            LiftOver_hal(
                Hal2chain_to_human.out.mix(Hal2chain_from_human.out).groupTuple(size:2),
                Human_bed
            )
            Evo_intersect( 
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
            Evo_intersect( 
                LiftOver.out[1].collect(),
                Human_bed,
                file("${baseDir}/bin/intersect_elements.py")
            )
            }
        }
        else if (params.CHAINS) { 
            Evo_match(
            file(params.tebag_db),
            file(params.quantification_file),
            file("${baseDir}/bin/tebag_match.py")
            )
            Upset_plot(
            file("${baseDir}/bin/create_upset_plot.R"),
            Evo_match.out
            )
        }
    if ( params.permutation_test ) { 
        bigbeds = Channel
        //.fromPath('s3://human-pangenomics/T2T/CHM13/assemblies/annotation/regulation/ENCODE/macs2_peak/*.bb')
        .fromPath("${params.bigbed_files}/*.bb")
            //.filter { filePath -> filePath.name.endsWith('.bb') }
            .map { filePath -> 
                def fileName = filePath.name
                def groupKey = fileName.split('\\.')[2] // Extract the value before .bb
                [groupKey, filePath]
            }.groupTuple()
    
            age_annotated_bed = Channel
            .fromPath(params.age_annotated_bed)
        
    splitbed_output = Splitbed(age_annotated_bed)
//splitbed_output.split_bed.flatten().combine(bigbeds).view()


    Runpermutation(
        bigbeds.combine(splitbed_output.split_bed.flatten()),
        //splitbed_output.split_bed.flatten().combine(bigbeds),
        file(params.fai)
    )


    }
}
#!/bin/bash

# lift back and forth between two organisms 
# args (in order): organism1 organism2 chain_directory bed_directory

organism_one=$1
organism_two=$2
chain_directory=$3
bed_directory=$4

echo $chain_directory/$organism_two-organism_one.chain.gz 
echo  $bed_directory/$organism_one.rmsk.bed 
#lift from organism one to organism two 
# organism_two.first.lifted.bed is the lifted coordinates to the new genome
# organism_two.first.unlifted.bed is anything that didn't lift 

liftOver $bed_directory/$organism_one.rmsk.bed \
	$chain_directory/$organism_two-$organism_one.chain.gz \
	$organism_one.$organism_two.first.lifted.bed \
	$organism_one.$organism_two.first.unlifted.bed

# lifts the new coordinates back to the first organism to check if lifting was successful 
liftOver $organism_one.$organism_two.first.lifted.bed \
	$chain_directory/$organism_one-$organism_two.chain.gz \
	$organism_one.$organism_two.check.lifted.bed \
	$organism_one.$organism_two.check.unlifted.bed

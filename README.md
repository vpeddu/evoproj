# TEBAG
## Transposable Elements By Age Group

## To Do's and Notes  

### 2022_07_07  

- Finish drafting liftover  
- Intersection code -- python or R?  
- Input: names, abundance, DESeq  
- VCF: think about  
- Test: LTR7up paper for comparison  
	- [link to test bed](https://github.com/LumpLord/Mosaic-cis-regulatory-evolution-drives-transcriptional-partitioning-of-HERVH-endogenous-retrovirus../blob/main/Aging/Liftover_aging/all_named.bed)
	- take this, use to filter our bed, use as test
- Paper: re-process some data, demonstrate performance  
- Undergrads: how to get involved

### sample run command 
```
nextflow run main.nf --generate_db --species_paths bin/sample_generate_csv.csv --human_bed ../beds/ucsc.rmsk.salmon.bed -with-docker ubuntu:18.04 -resume 
```
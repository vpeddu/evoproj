#!/bin/bash
    ls -lah

    # Process bed file with a simpler approach
    mkdir -p under_threshold
    
 awk '
BEGIN { FS="\t" }  # Set field separator to tab
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
}' test.bed

rm *\)n*

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

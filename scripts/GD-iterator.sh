#!/bin/bash

# Assign input arguments to variables
input_folder="$1"
starter="$2"
output_folder="$3"

starter_basename=$(basename "$starter" | sed 's/\.candidates\.fasta$//')

# Copy the starter file to the output file
cp "$starter" "$output_folder/GD-final.fa"
samtools faidx "$output_folder/GD-final.fa"

# Loop through all .fasta files in the input folder
for file in "$input_folder"/*.fasta; do
    # Get just the filename without the path
    filename=$(basename "$file" | sed 's/\.candidates\.fasta$//')
    echo "$filename"
    samtools faidx "$file"

    # Perform the BLAST operation
    blastn -query "$output_folder/GD-final.fa" -subject "$file" -outfmt 6 -out "${output_folder}/${filename}.blast.txt"

    # Check if the filename is not equal to the starter file
    if [ "$filename" != "$starter_basename" ]; then
        # Call Rscript (replace with actual script and arguments)
        # Rscript path_to_script.R arg1 arg2 ...
        Rscript "find_new_candidates.R" "${output_folder}/${filename}.blast.txt" "$output_folder/GD-final.fa.fai" "${file}.fai" "${output_folder}/${filename}.selected.txt"
        python "filter-fasta.py" "${file}" "${output_folder}/${filename}.selected.txt" "${file}.selected.fa"
        cat "${file}.selected.fa" >> "$output_folder/GD-final.fa"
        samtools faidx "$output_folder/GD-final.fa"
    fi
done

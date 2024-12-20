#!/bin/bash

# Parse command line arguments
while getopts "i:m:o:" opt; do
    case $opt in
        i) input_gtf="$OPTARG";;      # original GTF for chrM
        m) modified_gtf="$OPTARG";;   # modified GTF for non-chrM
        o) output_gtf="$OPTARG";;
        *) echo "Usage: $0 -i <input_gtf[.gz]> -m <modified_gtf> -o <output_gtf>" >&2
          exit 1;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_gtf" ] || [ -z "$modified_gtf" ] || [ -z "$output_gtf" ]; then
    echo "Usage: $0 -i <input_gtf[.gz]> -m <modified_gtf> -o <output_gtf>"
    exit 1
fi

# Check if input files exist
if [ ! -f "$input_gtf" ]; then
    echo "Input file $input_gtf does not exist!"
    exit 1
fi
if [ ! -f "$modified_gtf" ]; then
    echo "Modified GTF file $modified_gtf does not exist!"
    exit 1
fi

echo "Processing GTF file: $input_gtf"

# Create temporary directory
temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT

# Create header file
cat <<-'EOF' > "$temp_dir/header.txt"
#gtf-version 2.2
#!genome-build mCalJa1.2.pat.X
#!genome-build-accession NCBI_Assembly:GCF_011100555.1
#!annotation-date 03/02/2023
#!annotation-source NCBI RefSeq GCF_011100555.1-RS_2023_03
###
EOF

# Check if input is gzipped
if [[ "$input_gtf" == *.gz ]]; then
    cat_cmd="gunzip -c"
else
    cat_cmd="cat"
fi

# Process non-chrM entries from modified GTF
$cat_cmd "$modified_gtf" | grep -v "^chrM" > "$temp_dir/non_chrm.gtf"

# Process chrM entries from original GTF
$cat_cmd "$input_gtf" | \
    grep "^chrM" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {
        # Store the original attributes
        attr=$9
        
        # Remove quotes
        gsub(/"/, "", attr)
        
        # Process each attribute
        n=split(attr, attrs, ";")
        new_attrs=""
        gene_id=""
        has_gene_name=0
        
        for(i=1; i<=n; i++) {
            # Skip empty fields
            if(attrs[i] ~ /^[[:space:]]*$/) continue
            
            # Split into key-value
            split(attrs[i], kv, " ")
            key=kv[1]
            val=kv[2]
            
            # Process ID fields
            if(key ~ /_id$/) {
                if(index(val, ".") > 0) {
                    # Split version number
                    split(val, ver, ".")
                    val=ver[1]
                    ver_key=substr(key, 1, length(key)-2) "_version"
                    new_attrs = new_attrs "; " ver_key " \"" ver[2] "\""
                }
            }
            
            # Track if we have gene_name
            if(key == "gene_name") has_gene_name=1
            
            # Store gene value
            if(key == "gene") gene_val=val
            
            # Add to new attributes
            new_attrs = new_attrs "; " key " \"" val "\""
        }
        
        # Add gene_name if missing and we have gene
        if(!has_gene_name && gene_val != "") {
            new_attrs = new_attrs "; gene_name \"" gene_val "\""
        }
        
        # Remove leading separator
        sub(/^; /, "", new_attrs)
        
        # Output the modified line
        $9=new_attrs
        print
    }' > "$temp_dir/modified_chrm.gtf"

# Create MT gene annotations
awk -F'\t' '
$9 ~ /gene_id "MT-/ {
    split($9, attrs, ";")
    for (i in attrs) {
        if (attrs[i] ~ /gene_id/) {
            match(attrs[i], /"([^"]+)"/)
            gene_id = substr(attrs[i], RSTART+1, RLENGTH-2)
            
            if (!(gene_id in start)) {
                start[gene_id] = $4
                end[gene_id] = $5
                source[gene_id] = $2
                chrom[gene_id] = $1
                strand[gene_id] = $7
                attr[gene_id] = $9
            } else {
                if ($4 < start[gene_id]) start[gene_id] = $4
                if ($5 > end[gene_id]) end[gene_id] = $5
            }
        }
    }
}

END {
    for (gene_id in start) {
        print chrom[gene_id] "\t" \
              source[gene_id] "\t" \
              "gene" "\t" \
              start[gene_id] "\t" \
              end[gene_id] "\t" \
              "." "\t" \
              strand[gene_id] "\t" \
              "." "\t" \
              attr[gene_id]
    }
}' "$temp_dir/modified_chrm.gtf" > "$temp_dir/mt_gene_annotation.txt"

# Combine all parts into final output
cat "$temp_dir/header.txt" \
    "$temp_dir/non_chrm.gtf" \
    "$temp_dir/modified_chrm.gtf" \
    "$temp_dir/mt_gene_annotation.txt" > "$output_gtf"

echo "Done! Processed GTF file has been written to: $output_gtf" 
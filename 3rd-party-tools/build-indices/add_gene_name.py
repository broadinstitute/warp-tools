# Script for adding gene_name as an additional attribute
import re

input_gtf = "genes.gtf"
output_gtf = "output_gene_name.gtf"

with open(input_gtf, "r") as infile, open(output_gtf, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue

        fields = line.strip().split("\t")
        if len(fields) != 9:
            outfile.write(line)
            continue

        attributes = fields[8]

        # Extract gene value
        match = re.search(r'gene\s+"([^"]+)"', attributes)
        if match:
            gene_value = match.group(1)
            # Only add gene_name if it doesn't already exist
            if 'gene_name' not in attributes:
                attributes = attributes.replace(
                    f'gene "{gene_value}";',
                    f'gene "{gene_value}"; gene_name "{gene_value}";'
                )
        fields[8] = attributes
        outfile.write("\t".join(fields) + "\n")


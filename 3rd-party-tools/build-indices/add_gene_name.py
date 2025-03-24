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

        match = re.search(r'gene_id\s+"([^"]+)"', attributes)
        if match:
            gene_value = match.group(1)
            if 'gene_name' not in attributes:
                attributes = attributes.strip() + f' gene_name "{gene_value}";'

        fields[8] = attributes
        outfile.write("\t".join(fields) + "\n")

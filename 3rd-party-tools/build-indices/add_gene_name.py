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

        # Parse attributes into a dictionary
        attr_dict = dict(re.findall(r'(\S+)\s+"(.*?)"', attributes))

        # If gene_name is missing, try to add it
        if "gene_name" not in attr_dict:
            # Try to use 'gene' first, then 'gene_id', then 'transcript_id'
            gene_value = attr_dict.get("gene") or attr_dict.get("gene_id") or attr_dict.get("transcript_id")
            if gene_value:
                # Add the gene_name at the end of the attributes string
                attributes = attributes.strip() + f' gene_name "{gene_value}";'

        fields[8] = attributes
        outfile.write("\t".join(fields) + "\n")

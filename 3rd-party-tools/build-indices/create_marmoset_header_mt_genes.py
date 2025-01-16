import sys

def add_gene_entries(gtf_file):

    genes = {}
    
    # First collect ALL features for each gene
    print("Colleting features of gene")
    for line in open(gtf_file):
        fields = line.strip().split('\t')
        if len(fields) != 9:
            continue
            
        attrs = {}
        for attr in fields[8].split('; '):
            if attr:
                key, value = attr.split(' ', 1)
                attrs[key] = value.strip('"')
                
        if 'gene_id' in attrs and attrs['gene_id'].startswith('MT-'):
            gene_id = attrs['gene_id']
            # Track coordinates regardless of feature type
            if gene_id not in genes:
                genes[gene_id] = {
                    'chrom': fields[0],
                    'source': fields[1],
                    'min_start': int(fields[3]),
                    'max_end': int(fields[4]),
                    'strand': fields[6],
                    'attrs': attrs,
                    'features': set([fields[2]])  # Track feature types
                }
            else:
                genes[gene_id]['min_start'] = min(genes[gene_id]['min_start'], int(fields[3]))
                genes[gene_id]['max_end'] = max(genes[gene_id]['max_end'], int(fields[4]))
                genes[gene_id]['features'].add(fields[2])

    gene_entries = []
    for gene_id, info in sorted(genes.items()):
        attrs_str = f'gene_id "{gene_id}"; '
        if 'gene' in info['attrs']:
            attrs_str += f'gene "{info["attrs"]["gene"]}"; '
        if 'transcript_biotype' in info['attrs']:
            attrs_str += f'transcript_biotype "{info["attrs"]["transcript_biotype"]}"; '
        if 'gene_version' in info['attrs']:
            attrs_str += f'gene_version "{info["attrs"]["gene_version"]}"; '
        if 'gene_name' in info['attrs']:
            attrs_str += f'gene_name "{info["attrs"]["gene_name"]}"'
            
        gene_entry = [
            info['chrom'],
            info['source'],
            'gene',
            str(info['min_start']),
            str(info['max_end']),
            '.',
            info['strand'],
            '.',
            attrs_str
        ]
        if gene_entry[-1].endswith(';"; '):
            gene_entry[-1] = gene_entry[-1][:-3]
        gene_entries.append('\t'.join(gene_entry))
        


    
    print("Setting up header")
    with open(gtf_file) as f:

                # First print the new header
        new_header = """#gtf-version 2.2
#!genome-build mCalJa1.2.pat.X
#!genome-build-accession NCBI_Assembly:GCF_011100555.1
#!annotation-date 03/02/2023
#!annotation-source NCBI RefSeq GCF_011100555.1-RS_2023_03
##gff-version 2
##source-version rtracklayer 1.58.0
##date 2023-09-11"""
        print(new_header)

        for line in f:
            if not line.startswith("#"):
                print(line.strip())
    print("Writing gene entries")
    for entry in gene_entries:
        print(entry)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py input.gtf")
        sys.exit(1)
    add_gene_entries(sys.argv[1])



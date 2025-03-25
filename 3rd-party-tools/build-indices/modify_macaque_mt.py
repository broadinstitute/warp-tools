# Script to find mitochondrial genes and add them as gene entries
import sys

def load_mito_gene_list(txt_file):
    with open(txt_file) as f:
        return set(line.strip() for line in f if line.strip() and not line.startswith("#"))

def add_gene_entries(gtf_file, mito_txt_file):
    mito_genes = load_mito_gene_list(mito_txt_file)
    genes = {}

    # First collect features for mitochondrial genes only
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            attrs = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr:
                    key_value = attr.split(' ', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        attrs[key] = value.strip('"')

            gene_name = attrs.get("gene_name", "")
            gene_id = attrs.get("gene_id", "")
            if gene_name in mito_genes:
                if gene_id not in genes:
                    genes[gene_id] = {
                        'chrom': fields[0],
                        'source': fields[1],
                        'min_start': int(fields[3]),
                        'max_end': int(fields[4]),
                        'strand': fields[6],
                        'attrs': attrs,
                        'features': set([fields[2]])
                    }
                else:
                    genes[gene_id]['min_start'] = min(genes[gene_id]['min_start'], int(fields[3]))
                    genes[gene_id]['max_end'] = max(genes[gene_id]['max_end'], int(fields[4]))
                    genes[gene_id]['features'].add(fields[2])

    # Create synthetic 'gene' entries
    gene_entries = []
    for gene_id, info in sorted(genes.items()):
        attrs = info['attrs']
        attrs_str = f'gene_id "{gene_id}"; '
        if 'gene' in attrs:
            attrs_str += f'gene "{attrs["gene"]}"; '
        if 'transcript_biotype' in attrs:
            attrs_str += f'transcript_biotype "{attrs["transcript_biotype"]}"; '
        if 'gene_version' in attrs:
            attrs_str += f'gene_version "{attrs["gene_version"]}"; '
        if 'gene_name' in attrs:
            attrs_str += f'gene_name "{attrs["gene_name"]}";'

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
        gene_entries.append('\t'.join(gene_entry))

    # Output header and full GTF content + added gene entries
    with open(gtf_file) as f:
        for line in f:
            print(line.strip())
    for entry in gene_entries:
        print(entry)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 script.py input.gtf mito_gene_list.txt")
        sys.exit(1)
    add_gene_entries(sys.argv[1], sys.argv[2])


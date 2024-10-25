#!/usr/bin/env python3

import argparse
import gzip
import re

def get_ref_source(input_gtf):
    """Determine if the GTF is from RefSeq or GENCODE."""
    reference_source = ""
    with open(input_gtf, 'r') as fh:
        for line_num, line in enumerate(fh):
            if ('gene_type' in line) or ('transcript_type' in line):
                reference_source = "Gencode"
                print('Gencode Reference')
                break
            elif ('gene_biotype' in line) or ('transcript_biotype' in line):
                reference_source = "Refseq"
                print('Refseq reference')
                break
    return reference_source 

def is_valid_chromosome_region(fields, species):
    """
    Check if the feature is in a valid chromosome region.
    Only filters PAR regions for human.
    """
    if species != "human":
        return True
        
    chromosome = fields[0]
    if chromosome != "chrY":
        return True
        
    start_pos = int(fields[3])
    # For human chrY, check if it's in the valid region (2752083-56887903)
    return 2752083 <= start_pos < 56887903

def is_allowable_biotype(biotype, ref_source, features_dic=None):
    """Check if biotype should be included based on reference source."""
    # Handle mitochondrial genes first
    if features_dic and (features_dic.get('gene', '').startswith('MT-') or 
                        features_dic.get('gene_name', '').startswith('MT-')):
        return True

    if ref_source == "Refseq":
        allowed_types = {
            'mRNA',           
            'tRNA',
            'rRNA',
            'C_gene_segment',
            'V_gene_segment',
            'lnc_RNA',      
        }
    else:  # Gencode/Ensembl
        allowed_types = {
            'protein_coding',
            'protein_coding_LoF',
            'lncRNA',
            'tRNA',
            'rRNA',
            'mRNA',
            'IG_C_gene',
            'IG_D_gene',
            'IG_J_gene',
            'IG_LV_gene',
            'IG_V_gene',
            'IG_V_pseudogene',
            'IG_J_pseudogene',
            'IG_C_pseudogene',
            'TR_C_gene',
            'TR_D_gene',
            'TR_J_gene',
            'TR_V_gene',
            'TR_V_pseudogene',
            'TR_J_pseudogene'
        }
    
    # Handle cases where biotype is missing but it's an MT gene
    if biotype is None and features_dic and features_dic.get('gene', '').startswith('MT-'):
        return True
        
    return biotype in allowed_types

def get_features(features):
    """Parse GTF attributes into a dictionary."""
    features_dict = {}
    key_value_pairs = features.split(';')
    for pair in key_value_pairs:
        key_value = pair.strip().split(' ', 1)
        if len(key_value) == 2:
            key, value = key_value
            key = key.strip()
            value = value.strip(' "')
            if key:
                features_dict[key] = value
        elif len(key_value) == 1:
            key = key_value[0].strip()
            features_dict[key] = ""
    return features_dict

def modify_attr(features_dic):
    """Modify attribute fields to version format."""
    for key in features_dic.copy():
        if key in ["exon_id", "gene_id", "transcript_id"]:
            version_key = key.replace("_id", "_version")
            if '.' in features_dic[key]:
                features_dic[key], version = features_dic[key].split(".", 1)
            else:
                features_dic[key] = features_dic[key]
                version = 0

            if version_key not in features_dic:
                features_dic[version_key] = version
        if key == 'gene' and "gene_name" not in features_dic.keys():
            features_dic['gene_name'] = features_dic['gene']

    modified_features = "; ".join([f'{key} "{value}"' for key, value in features_dic.items() if key and value is not None])
    return modified_features

def get_gene_ids_Refseq(input_gtf, species):
    """Collect gene IDs from RefSeq GTF that meet filtering criteria."""
    gene_ids = []
    with open(input_gtf, 'r') as input_file:
        for line in input_file:
            if not line.startswith('#'):
                fields = [x.strip() for x in line.strip().split('\t')]
                if not is_valid_chromosome_region(fields, species):
                    continue
                    
                if fields[2] == 'transcript':
                    features = re.sub('"', '', line.strip().split('\t')[8].strip())
                    features_dic = get_features(features)

                    # Special handling for mitochondrial genes
                    if fields[0] == "chrM" or features_dic.get('gene', '').startswith('MT-'):
                        gene = features_dic['gene_id'].split('.', 1)[0]
                        if gene not in gene_ids:
                            gene_ids.append(gene)
                        continue

                    biotype = features_dic.get('transcript_biotype') or features_dic.get('gene_biotype')
                    if biotype and is_allowable_biotype(biotype, "Refseq", features_dic):
                        gene = features_dic['gene_id'].split('.', 1)[0]
                        if gene not in gene_ids:
                            gene_ids.append(gene)
                    elif features_dic.get('gene', '').startswith('MT-'):
                        # Handle MT genes without biotype
                        gene = features_dic['gene_id'].split('.', 1)[0]
                        if gene not in gene_ids:
                            gene_ids.append(gene)
    return gene_ids

def get_gene_ids_Gencode(input_gtf, species):
    """Collect gene IDs from GENCODE GTF that meet filtering criteria."""
    gene_ids = set()
    with open(input_gtf, 'r') as input_file:
        print("open successful")
        for line in input_file:
            if line.startswith('#'):
                continue
            fields = [x.strip() for x in line.strip().split('\t')]
            
            if not is_valid_chromosome_region(fields, species):
                continue
                
            if fields[2] != 'transcript':
                continue
            
            features = re.sub('"', '', line.strip().split('\t')[8].strip())
            features_dic = get_features(features)

            # Special handling for mitochondrial genes
            if fields[0] == "chrM" or features_dic.get('gene', '').startswith('MT-'):
                gene = features_dic['gene_id']
                if gene not in gene_ids:
                    gene_ids.add(gene)
                continue

            gene_type = features_dic.get('gene_type')
            transcript_type = features_dic.get('transcript_type') or features_dic.get('transcript_biotype')
            
            # Handle MT genes without biotype
            if features_dic.get('gene', '').startswith('MT-'):
                gene = features_dic['gene_id']
                if gene not in gene_ids:
                    gene_ids.add(gene)
                continue
                
            if (gene_type and is_allowable_biotype(gene_type, "Gencode", features_dic) and
                (transcript_type is None or is_allowable_biotype(transcript_type, "Gencode", features_dic))):
                gene = features_dic['gene_id']
                if gene not in gene_ids:
                    gene_ids.add(gene)

    return list(gene_ids)

def filter_gtf(input_gtf, output_gtf, species, gene_ids):
    """Filter and write the GTF file."""
    with open(input_gtf, 'r') as input_file:
        with open(output_gtf, 'w') as output_gtf:
            for line in input_file:
                if line.startswith("#"):
                    output_gtf.write(line)
                    continue
                    
                fields = [x.strip() for x in line.strip().split("\t")]
                
                if not is_valid_chromosome_region(fields, species):
                    continue
                    
                features = re.sub('"', '', line.strip().split('\t')[8].strip())
                if features.endswith(';'):
                    features = features[:-1]                     
                features_dic = get_features(features)
                
                # For human, skip XGY2 gene explicitly
                if species == "human" and features_dic.get('gene_id') == "ENSG00000290840":
                    continue
                    
                if features_dic['gene_id'] in gene_ids:
                    modified_fields = fields.copy()
                    modified_fields[8] = modify_attr(features_dic)
                    output_gtf.write("{}".format("\t".join(modified_fields)+ "\n"))

def main():
    """Main function to process GTF files."""
    parser = argparse.ArgumentParser(description="Filter GTF files based on biotypes and format specifications")
    parser.add_argument(
        "--input-gtf",
        "-i",
        dest="input_gtf",
        required=True,
        help="input GTF file"
    )
    parser.add_argument(
        "--output-gtf",
        "-o",
        dest="output_gtf",
        required=True,
        help="output GTF file"
    )
    parser.add_argument(
        "--species",
        "-s",
        dest="species",
        required=True,
        help="Species name (e.g., human, mouse, etc.)"
    )
    args = parser.parse_args()

    print(f"Processing {args.species} GTF file...")

    print("Determining reference source...")
    ref_source = get_ref_source(input_gtf=args.input_gtf)
    
    print("Collecting gene IDs...")
    if ref_source == "Refseq":
        gene_ids = get_gene_ids_Refseq(input_gtf=args.input_gtf, species=args.species)
    elif ref_source == "Gencode":
        gene_ids = get_gene_ids_Gencode(input_gtf=args.input_gtf, species=args.species)
    else:
        print("Error: The input GTF file format is not recognized. Must be either GENCODE or RefSeq format.")
        exit(1)

    print(f"Found {len(gene_ids)} genes matching criteria")
    print("Writing filtered GTF file...")
    
    filter_gtf(args.input_gtf, args.output_gtf, args.species, gene_ids)
    print("Done! Filtered GTF file has been written to:", args.output_gtf)

if __name__ == "__main__":
    main()
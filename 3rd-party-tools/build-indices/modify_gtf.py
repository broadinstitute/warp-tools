import argparse
import gzip
import re

def get_ref_source(input_gtf):
    reference_source = ""
    with open(input_gtf, 'r') as fh:
        for line_num, line in enumerate(fh):
            # search string
            if ('gene_type' in line) or ('transcript_type' in line):
                reference_source = "Gencode"
                print('Gencode Reference')
                # don't continue the search if found at least one instance
                break
            elif ('gene_biotype' in line) or ('transcript_biotype' in line):
                reference_source = "Refseq"
                print('Refseq reference')
                # don't continue the search if found at least one instance
                break
    return reference_source 

def get_biotypes(biotypes_file_path):
    allowable_biotypes= []
    with open(biotypes_file_path,'r',encoding='utf-8-sig') as biotypesFile:
        for line in biotypesFile:
            if not line.startswith("#"):
                fields = [x.strip() for x in line.strip().split("\t")]
                if fields[1] == "Y" or fields[1]=="y":
                    allowable_biotypes.append(fields[0])
    return  allowable_biotypes

def get_features(features):
    features_dict = {}
    key_value_pairs = features.split(';')
    #Process each key-value pair
    for pair in key_value_pairs:
        # Split each pair into key and value
        key_value = pair.strip().split(' ', 1)

        # Ensure there is a key (and at least one element in the key-value pair)
        if len(key_value) == 2:
            key, value = key_value
            key = key.strip()
            value = value.strip(' "')

            # Add the key-value pair to the dictionary
            if key:
                features_dict[key] = value
        elif len(key_value) == 1:
            # Handle the case where a key is present but the value is missing
            key = key_value[0].strip()
            features_dict[key] = ""
    return features_dict

def modify_attr(features_dic):
    modified_features = ""
    if "gene_name" is not in features_dic.keys():
        if "gene" in features_dic.keys():
            features_dic["gene_name"] = features_dic["gene"]
        else:
            features_dic["gene_name"] = features_dic["gene_id"]

    for key in features_dic.copy():
        if key in ["exon_id","gene_id","transcript_id"]:
            version_key = key.replace("_id", "_version")
            # Check if the gene_id has a version and split on the period
            if '.' in features_dic[key]:
                features_dic[key], version = features_dic[key].split(".",1)
            else:
                features_dic[key] = features_dic[key]
                version = 0

            # If the version ids are not present in the dictionary, add them
            if version_key not in features_dic:
                features_dic[version_key] = version

    modified_features = "; ".join([f'{key} "{value}"' for key, value in features_dic.items() if key and value is not None])

    return modified_features

def get_gene_ids_Refseq(input_gtf, biotypes):
    gene_ids =[]
    with open(input_gtf, 'r') as input_file:
        for line in input_file:
            if not line.startswith('#'):
                fields = [x.strip() for x in line.strip().split('\t')]
                if fields[2]=='gene':
                    features = re.sub('"', '', line.strip().split('\t')[8].strip())
                    features_dic = get_features(features)
                    if ('gene_biotype' in features_dic.keys()) and (features_dic['gene_biotype'] in biotypes):
                        gene=features_dic['gene_id'].split('.',1)[0]
                        if gene not in gene_ids:
                            gene_ids.append(gene)
    input_file.close()
    return gene_ids

def get_gene_ids_Gencode(input_gtf, biotypes):
    gene_ids = set()
    with open(input_gtf, 'r') as input_file:
        print("open succesful")
        for line in input_file:
            if line.startswith('#'):
                continue
            fields = [x.strip() for x in line.strip().split('\t')]
            if fields[2] != 'transcript':
                continue
            features = re.sub('"', '', line.strip().split('\t')[8].strip())
            features_dic = get_features(features)

            if (features_dic['gene_type'] not in biotypes) or (features_dic['transcript_type'] not in biotypes):
                continue

            if 'tag' in features_dic:
                if ('readthrough_transcript' not in features_dic['tag']) and (
                    'PAR' not in features_dic['tag']):
                    gene=features_dic['gene_id']
                    if gene not in gene_ids:
                        gene_ids.add(gene)

            else:
                gene=features_dic['gene_id']
                if gene not in gene_ids:
                    gene_ids.add(gene)

    input_file.close()
    return list(gene_ids)

def main():
    """ This script filters a GTF file for all genes that have at least one transcript with a biotype.
    The complete list of biotypes is passed and the desired biotypes are marked by a boolean Y or N.
    The new gtf file attributes exon_id, gene_id and transcript_id are modified by removing the versions to match Cellranger IDs.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-gtf",
        "-i",
        dest="input_gtf",
        default=None,
        required=True,
        help="input GTF",
    )
    parser.add_argument(
        "--output-gtf",
        "-o",
        dest="output_gtf",
        default=None,
        help="output GTF"
    )
    parser.add_argument(
        "--biotypes",
        "-b",
        dest="biotypes_file",
        default=None,
        required=True,
        help="List of all gene_type or transcript_type fields and a boolean to choose which biotypes to filter in the GTF file."
    )
    args = parser.parse_args()

    allowable_biotypes = get_biotypes(biotypes_file_path=args.biotypes_file)
    print("Filtering the gtf for biotypes:", allowable_biotypes)
    ref_source = get_ref_source(input_gtf=args.input_gtf)
    if ref_source == "Refseq":
        gene_ids = get_gene_ids_Refseq(input_gtf=args.input_gtf, biotypes=allowable_biotypes)
    elif ref_source == "Gencode":
        gene_ids = get_gene_ids_Gencode(input_gtf=args.input_gtf, biotypes=allowable_biotypes)
    else:
        print("The input gtf file cannot be processed with this script. Provided biotype attributes should be from Gencode or Refseq.")
        exit(0)

    print("Writing the output file based on selected genes with biotype")
    with open(args.input_gtf, 'r') as input_file:
        with open(args.output_gtf, 'w') as output_gtf:
            for line in input_file:
                if line.startswith("#"):
                    output_gtf.write(line.strip() + "\n")
                else:
                    fields = [x.strip() for x in line.strip().split("\t")]
                    features = re.sub('"', '', line.strip().split('\t')[8].strip())
                    # Remove the trailing semicolon if it exists
                    if features.endswith(';'):
                      features = features[:-1]                     
                    #print("Features:",features,'\n', fields[8])
                    features_dic = get_features(features)
                    #print(features_dic)
                    if features_dic['gene_id'] in gene_ids:
                        # The two lines below for modified filed were moved into this if statement
                        # We want to find the valid genes first and then modify the GTF fields
                        modified_fields = fields.copy()
                        modified_fields[8] = modify_attr(features_dic)
                     #   print("modfied:",modified_fields[8])
                        output_gtf.write("{}".format("\t".join(modified_fields)+ "\n"))

if __name__ == "__main__":
    main()

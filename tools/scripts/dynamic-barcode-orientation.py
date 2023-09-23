import argparse
def read_whitelist(whitelist_file):
    # Read and return the whitelist from a file (one barcode per line)
    with open(whitelist_file, 'r') as file:
        return set(line.strip() for line in file)

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def count_matching_barcodes(input_file, whitelist_file):
    # Read the whitelist into a set
    whitelist = read_whitelist(whitelist_file)
    
    # Initialize variables to count matched barcodes
    count_forward_start = 0
    count_reverse_comp_start = 0
    count_forward_end = 0
    count_reverse_comp_end = 0
    
    with open(input_file, 'r') as file:
        # Read the barcodes from the input file
        barcodes = [line.strip() for line in file.readlines()]
        
        for barcode in barcodes:
            # Get reverse complement of entire barcode sequence
            reverseComplement = reverse_complement(barcode)
            
            # Match the first 16 base pairs in the forward direction
            start_forward = barcode[:16]
            if start_forward in whitelist:
                count_forward_start += 1
            
            # Match the first 16 base pairs in the reverse complement direction
            start_reverse_comp = reverseComplement[:16]
            if start_reverse_comp in whitelist:
                count_reverse_comp_start += 1
            
            # Match the last 16 base pairs in the forward direction
            end_forward = barcode[-16:]
            if end_forward in whitelist:
                count_forward_end += 1
            
            # Match the last 16 base pairs in the reverse complement direction
            end_reverse_comp = reverseComplement[-16:]
            if end_reverse_comp in whitelist:
                count_reverse_comp_end += 1
    
    return (
        count_forward_start,
        count_reverse_comp_start,
        count_forward_end,
        count_reverse_comp_end
    )

def determine_best_matching_method(
    count_forward_start,
    count_reverse_comp_start,
    count_forward_end,
    count_reverse_comp_end
):
    counts = [
        ("FIRST_BP", count_forward_start),
        ("FIRST_BP_RC", count_reverse_comp_start),
        ("LAST_BP", count_forward_end),
        ("LAST_BP_RC", count_reverse_comp_end)
    ]
    
    # Find the matching method with the highest count
    best_match, highest_count = max(counts, key=lambda x: x[1])
    if highest_count < 50:
        raise Exception("Less than one percent of barcodes match whitelist")
    return best_match

def main():
    parser = argparse.ArgumentParser(description="Count matching DNA barcodes and determine the best matching method.")
    parser.add_argument("input_file", help="Path to the input file containing DNA barcodes.")
    parser.add_argument("whitelist_file", help="Path to the whitelist file containing allowed barcodes.")
    parser.add_argument("output_file", help="Path to the output file for the best matching method parameter.")

    args = parser.parse_args()
    # Define input and output files here you remove argparse
    #input_file = "downsample.fq"  # Change to the path of your input file
    #whitelist_file = "~{whitelist}"  # Change to the path of your whitelist file
    #output_file = "best_match.txt"  # Change to the path of your output file

    (
        count_forward_start,
        count_reverse_comp_start,
        count_forward_end,
        count_reverse_comp_end
    ) = count_matching_barcodes(args.input_file, args.whitelist_file)

    print("Number of Forward Matches (Start 16 bp):", count_forward_start)
    print("Number of Reverse Complement Matches (Start 16 bp):", count_reverse_comp_start)
    print("Number of Forward Matches (End 16 bp):", count_forward_end)
    print("Number of Reverse Complement Matches (End 16 bp):", count_reverse_comp_end)

    best_matching_method = determine_best_matching_method(
        count_forward_start,
        count_reverse_comp_start,
        count_forward_end,
        count_reverse_comp_end
    )
    print(best_matching_method)

    with open(args.output_file, 'w') as outfile:
        outfile.write(best_matching_method)

if __name__ == "__main__":
    main()        

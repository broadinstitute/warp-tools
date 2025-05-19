import anndata
import sys

def load_anndata(file_path):
    # Load an AnnData object from the specified H5AD file
    return anndata.read_h5ad(file_path)

def main():
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <test_h5ad_file> <truth_h5ad_file>")
        sys.exit(1)

    # Get the input file paths from command-line arguments
    test_h5ad_path = sys.argv[1]
    truth_h5ad_path = sys.argv[2]

    # Load AnnData objects from the H5AD files
    test = load_anndata(test_h5ad_path)
    truth = load_anndata(truth_h5ad_path)

    # Print some information about the loaded AnnData objects (optional)
    print("Loaded test AnnData:")
    print(test)

    print("\nLoaded truth AnnData:")
    print(truth)

    # Now you can use the 'test' and 'truth' AnnData objects as needed in your script
    # Do comparisons between test and truth obs
    print("Comparing test cell metrics to truth cell metrics using truth as ref")
    command1="test.obs."
    command2="truth.obs."
    for x in truth.obs.columns:
        z = eval(command1+x)
        y = eval(command2+x)
        if z.equals(y)==False:
            print("Cell Metric Column does not match:")
            print(x)
    print("Comparing test gene metrics to truth gene metrics using truth as ref")
    command3="test.var."
    command4="truth.var."
    for x in truth.var.columns:
        z = eval(command3+x)
        y = eval(command4+x)
        if z.equals(y)==False:
            print("Gene Metric Column does not match:")
            print(x)
    print("Testing matrix count sums")
    if test.X.sum()==truth.X.sum():
        print("Counts match")
    else:
        print("Counts do not match")
    print("Testing that multimapped_reads + unique_reads = n_reads in test")
    if test.obs.n_reads.sum()==(test.obs.reads_mapped_multiple.sum()+test.obs.reads_mapped_uniquely.sum()):
        print("N_reads = reads_mapped_multiple+unique")
    else:
        print("Sums of unique, multimapped do not equal n_reads")

if __name__ == "__main__":
    main()

import pytest
import pandas as pd
from ..scripts.combine_shard_metrics import merge_matrices

# Mock data setup for summary file
@pytest.fixture
def summary_file(tmpdir):
    data = {'metric': ['Number of Reads', 'Another Metric'], 'value': [100, 200], 'shard': ['shard1', 'shard1']}
    df = pd.DataFrame(data)
    file = tmpdir.join("./Human_v2_snRNA_summary.txt")
    df.to_csv(file, index=False, header=None)
    return str(file)

# Mock data setup for align file
@pytest.fixture
def align_file(tmpdir):
    data = {'metric': ['Reads Mapped to Genome: Unique', 'Q30 Bases in RNA read'], 'value': [50, 75], 'shard': ['shard1', 'shard1']}
    df = pd.DataFrame(data)
    file = tmpdir.join("./input/Human_v2_snRNA_align_features.txt")
    df.to_csv(file, sep="\t", index=False, header=None)  # Assuming TSV format
    return str(file)

# Mock data setup for cell reads file
@pytest.fixture
def cell_reads_file(tmpdir):
    data = {'CB': ['CB1', 'CB2'], 'exonicAS': [10, 20], 'intronicAS': [5, 15], 'exonic': [20, 30]}
    df = pd.DataFrame(data)
    file = tmpdir.join("/input/Human_v2_snRNA_cell_reads.tsv")
    df.to_csv(file, sep="\t", index=False, header=None)
    return str(file)

# Test function for merge_matrices
def test_merge_matrices(summary_file, align_file, cell_reads_file):
    counting_mode = "sc_rna"
    uniform_barcodes = "./input/barcodes.tsv"
    uniform_mtx = "./input/matrix.mtx"
    df = merge_matrices(summary_file, align_file, cell_reads_file, counting_mode, uniform_barcodes, uniform_mtx)
    
    # Asserting expected dataframe size
    assert df.shape[0] == 1  # Expecting one row of results
    assert df.shape[1] > 1  # Expecting multiple columns of metrics
    
    # Asserting expected values for specific metrics
    # Adjust the keys according to your actual expected output
    assert df['number_of_reads'][0] == 100
    assert df['estimated_cells'][0] == 2  # Assuming 'estimated_cells' is a key and the expected result
    # Add more assertions here based on expected results
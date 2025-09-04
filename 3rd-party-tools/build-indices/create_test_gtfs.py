# test_data/test1.gtf
import os

def create_test_gtfs():
    # Create test_data directory if it doesn't exist
    os.makedirs("test_data", exist_ok=True)
    
    # Test GTF 1
    test_gtf1 = """#!genome-build GRCh38.p13
#!genome-version GRCh38
chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";
chr1\tENSEMBL\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";
chrM\tHAVANA\tgene\t577\t647\t.\t+\t.\tgene_id "MT-TF-1"; gene_name "MT-TF";"""

    # Test GTF 2 (with some differences)
    test_gtf2 = """#!genome-build GRCh38.p13
#!genome-version GRCh38
chr1\tHAVANA\tgene\t11869\t14410\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1";
chr1\tHAVANA\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";
chr1\tENSEMBL\texon\t11869\t12227\t.\t+\t.\tgene_id "ENSG00000223972"; transcript_id "ENST00000456328";
chrM\tHAVANA\tgene\t577\t647\t.\t+\t.\tgene_id "MT-TF-2"; gene_name "MT-TF-modified";"""

    # Write test files
    with open("test_data/test1.gtf", "w") as f:
        f.write(test_gtf1)
    
    with open("test_data/test2.gtf", "w") as f:
        f.write(test_gtf2)

if __name__ == "__main__":
    create_test_gtfs()
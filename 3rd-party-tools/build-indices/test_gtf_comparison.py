import unittest
import os
from compare_gtfs import compare_gtfs
import shutil

class TestGTFComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Create test GTFs before running tests
        from create_test_gtfs import create_test_gtfs
        create_test_gtfs()
        
        # Create output directory
        os.makedirs("test_output", exist_ok=True)

    def setUp(self):
        # Clean output directory before each test
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")
        os.makedirs("test_output")

    def test_gtf_comparison(self):
        # Run comparison
        compare_gtfs(
            "test_data/test1.gtf",
            "test_data/test2.gtf",
            "test_output/comparison"
        )
        
        # Check if output files were created
        self.assertTrue(os.path.exists("test_output/comparison_structural_diff.txt"))
        self.assertTrue(os.path.exists("test_output/comparison_attribute_diff.txt"))
        self.assertTrue(os.path.exists("test_output/comparison_gene_diff.txt"))
        
        # Read and check content of structural diff file
        with open("test_output/comparison_structural_diff.txt", "r") as f:
            content = f.read()
            self.assertIn("Total rows in GTF1: 4", content)
            self.assertIn("Total rows in GTF2: 4", content)
            
        # Read and check gene diff file
        with open("test_output/comparison_gene_diff.txt", "r") as f:
            content = f.read()
            self.assertIn("Total genes in GTF1: 2", content)
            self.assertIn("Total genes in GTF2: 2", content)
            self.assertIn("MT-TF-1", content)  # Check for gene in GTF1
            self.assertIn("MT-TF-2", content)  # Check for gene in GTF2
            self.assertIn("Genes only in GTF1: 1", content)
            self.assertIn("Genes only in GTF2: 1", content)

    @classmethod
    def tearDownClass(cls):
        # Clean up test data after all tests
        if os.path.exists("test_data"):
            shutil.rmtree("test_data")
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")

if __name__ == '__main__':
    unittest.main()
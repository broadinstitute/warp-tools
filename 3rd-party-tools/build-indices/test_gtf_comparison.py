import unittest
import os
import shutil
from compare_gtfs import compare_gtfs
from modify_gtf import modify_gtf  # Import the modify_gtf function if available

class TestGTFComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Create output directory
        os.makedirs("test_output", exist_ok=True)

    def setUp(self):
        # Clean output directory before each test
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")
        os.makedirs("test_output")

    def test_gtf_modification_and_comparison(self):
        input_gtf = "test_data/test1.gtf"
        modified_gtf = "test_output/modified_test1.gtf"
        biotypes_file = "Biotypes.tsv"
        
        # Verify input files exist
        self.assertTrue(os.path.exists(input_gtf), f"Input GTF file not found: {input_gtf}")
        self.assertTrue(os.path.exists(biotypes_file), f"Biotypes file not found: {biotypes_file}")
        
        # Run modification script (if modify_gtf is imported as a function)
        # If it's not imported as a function, you might need to use subprocess to run it
        import subprocess
        result = subprocess.run([
            'python', 'modify_gtf.py',
            '-i', input_gtf,
            '-o', modified_gtf,
            '-b', biotypes_file
        ], capture_output=True, text=True)
        
        self.assertEqual(result.returncode, 0, f"modify_gtf.py failed: {result.stderr}")
        
        # Verify modified file was created
        self.assertTrue(os.path.exists(modified_gtf))
        
        # Run comparison
        compare_gtfs(
            input_gtf,
            modified_gtf,
            "test_output/comparison"
        )
        
        # Check if output files were created
        self.assertTrue(os.path.exists("test_output/comparison_structural_diff.txt"))
        self.assertTrue(os.path.exists("test_output/comparison_attribute_diff.txt"))
        self.assertTrue(os.path.exists("test_output/comparison_gene_diff.txt"))
        
        # Read and check content of files
        with open("test_output/comparison_structural_diff.txt", "r") as f:
            content = f.read()
            self.assertIn("Total rows in GTF1", content)
            self.assertIn("Total rows in GTF2", content)
            
        with open("test_output/comparison_gene_diff.txt", "r") as f:
            content = f.read()
            self.assertIn("Total genes in GTF1", content)
            self.assertIn("Total genes in GTF2", content)

    @classmethod
    def tearDownClass(cls):
        # Clean up test output after all tests
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")

if __name__ == '__main__':
    unittest.main()
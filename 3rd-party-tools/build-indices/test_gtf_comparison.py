import unittest
import os
import sys
import shutil
import subprocess
from compare_gtfs import compare_gtfs, parse_attributes

class TestGTFModification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test directories and reference data"""
        # Create all required directories
        for directory in ["test_data", "test_output", "test_output/comparison_files", "test_data/reference_outputs"]:
            os.makedirs(directory, exist_ok=True)

    def setUp(self):
        """Clean and recreate test output directories"""
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")
        os.makedirs("test_output")
        os.makedirs("test_output/comparison_files")

    def run_modify_gtf(self, input_gtf, output_gtf, biotypes_file):
        """Run modify_gtf.py with given parameters"""
        result = subprocess.run([
            sys.executable,
            'modify_gtf.py',
            '-i', input_gtf,
            '-o', output_gtf,
            '-b', biotypes_file
        ], capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error output: {result.stderr}")
            print(f"Standard output: {result.stdout}")
        
        return result

    def compare_gtf_contents(self, file1, file2):
        """Compare GTF files with normalized content"""
        with open(file1) as f1, open(file2) as f2:
            lines1 = set(self.normalize_gtf_line(line) for line in f1)
            lines2 = set(self.normalize_gtf_line(line) for line in f2)
        
        only_in_1 = lines1 - lines2
        only_in_2 = lines2 - lines1
        
        return only_in_1, only_in_2

    def normalize_gtf_line(self, line):
        """Normalize a GTF line for comparison"""
        if line.startswith('#'):
            return line.strip()
        
        fields = line.strip().split('\t')
        if len(fields) != 9:
            return line.strip()
        
        # Parse and sort attributes to ensure consistent ordering
        attrs = parse_attributes(fields[8])
        sorted_attrs = '; '.join(f'{k} "{v}"' for k, v in sorted(attrs.items()))
        fields[8] = sorted_attrs
        
        return '\t'.join(fields)

    def create_summary_file(self, test_input, new_output, comparison_results):
        """Create a summary file of the test results"""
        summary_path = "test_output/comparison_files/test_summary.txt"
        with open(summary_path, 'w') as f:
            f.write("GTF Modification Test Summary\n")
            f.write("=========================\n\n")
            f.write(f"Input GTF: {test_input}\n")
            f.write(f"Output GTF: {new_output}\n\n")
            f.write("Comparison Results:\n")
            f.write("------------------\n")
            f.write(comparison_results)

    def test_gtf_modification_against_reference(self):
        """Test GTF modification against reference outputs"""
        test_input = "test_data/test1.gtf"
        reference_output = "test_data/reference_outputs/reference_modified.gtf"
        new_output = "test_output/new_modified.gtf"
        biotypes_file = "Biotypes.tsv"
        
        # Verify input files exist
        self.assertTrue(os.path.exists(test_input), f"Test input GTF not found: {test_input}")
        self.assertTrue(os.path.exists(biotypes_file), f"Biotypes file not found: {biotypes_file}")
        
        # Run current version of modify_gtf
        result = self.run_modify_gtf(test_input, new_output, biotypes_file)
        print("Script output:", result.stdout)
        print("Script errors:", result.stderr)
        
        self.assertEqual(result.returncode, 0, 
                        f"modify_gtf.py failed with error: {result.stderr}")
        
        # Verify the modified file was created
        self.assertTrue(os.path.exists(new_output), 
                       f"Modified GTF file was not created: {new_output}")
        
        # If reference output doesn't exist, create it
        if not os.path.exists(reference_output):
            print("Creating reference output for the first time...")
            os.makedirs(os.path.dirname(reference_output), exist_ok=True)
            shutil.copy(new_output, reference_output)
            self.skipTest("Reference output created. Run tests again to compare.")
        
        # Now run the comparison
        compare_gtfs(
            test_input,
            new_output,
            "test_output/comparison_files/comparison"
        )
        
        # Verify comparison output files were created
        expected_files = [
            "test_output/comparison_files/comparison_structural_diff.txt",
            "test_output/comparison_files/comparison_gene_diff.txt",
            "test_output/comparison_files/comparison_attribute_diff.txt"
        ]
        for file in expected_files:
            self.assertTrue(os.path.exists(file), f"Expected output file not created: {file}")
        
        # Compare normalized contents
        only_in_ref, only_in_new = self.compare_gtf_contents(reference_output, new_output)
        
        # Create detailed report regardless of differences
        report = ["GTF Comparison Results:"]
        if only_in_ref or only_in_new:
            if only_in_ref:
                report.append("\nLines only in reference GTF:")
                for line in sorted(only_in_ref)[:5]:
                    report.append(f"REF: {line}")
            
            if only_in_new:
                report.append("\nLines only in new GTF:")
                for line in sorted(only_in_new)[:5]:
                    report.append(f"NEW: {line}")
        else:
            report.append("\nNo differences found between reference and new GTF files.")
        
        # Always write the report
        with open("test_output/comparison_files/difference_report.txt", 'w') as f:
            f.write("\n".join(report))
        
        # Create summary file
        self.create_summary_file(test_input, new_output, "\n".join(report))
        
        # If there are differences, fail the test
        if only_in_ref or only_in_new:
            self.fail("\n".join(report[:10]) + 
                     "\n...\nSee test_output/comparison_files/difference_report.txt for full details")

    @classmethod
    def tearDownClass(cls):
        """Preserve outputs for GitHub Actions"""
        # Don't clean up - let GitHub Actions collect the artifacts
        pass

if __name__ == '__main__':
    unittest.main()
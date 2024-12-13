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
        os.makedirs("test_data", exist_ok=True)
        os.makedirs("test_output", exist_ok=True)
        os.makedirs("test_data/reference_outputs", exist_ok=True)
        
        if not os.path.exists("test_data/test1.gtf"):
            from create_test_gtfs import create_test_gtfs
            create_test_gtfs()

    def setUp(self):
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")
        os.makedirs("test_output")

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

    def compare_gtf_contents(self, file1, file2):
        """Compare GTF files with normalized content"""
        with open(file1) as f1, open(file2) as f2:
            lines1 = set(self.normalize_gtf_line(line) for line in f1)
            lines2 = set(self.normalize_gtf_line(line) for line in f2)
        
        only_in_1 = lines1 - lines2
        only_in_2 = lines2 - lines1
        
        return only_in_1, only_in_2

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
        
        # If reference output doesn't exist, create it
        if not os.path.exists(reference_output):
            print("Creating reference output for the first time...")
            os.makedirs(os.path.dirname(reference_output), exist_ok=True)
            shutil.copy(new_output, reference_output)
            self.skipTest("Reference output created. Run tests again to compare.")

        # Compare normalized contents
        only_in_ref, only_in_new = self.compare_gtf_contents(reference_output, new_output)
        
        # If there are real differences, create a detailed report
        if only_in_ref or only_in_new:
            report = ["Significant differences found between reference and new GTF:"]
            
            if only_in_ref:
                report.append("\nLines only in reference GTF:")
                for line in sorted(only_in_ref)[:5]:
                    report.append(f"REF: {line}")
            
            if only_in_new:
                report.append("\nLines only in new GTF:")
                for line in sorted(only_in_new)[:5]:
                    report.append(f"NEW: {line}")
            
            # Write detailed report
            with open("test_output/difference_report.txt", 'w') as f:
                f.write("\n".join(report))
            
            # Show beginning of the report in the failure message
            self.fail("\n".join(report[:10]) + 
                     "\n...\nSee test_output/difference_report.txt for full details")

    @classmethod
    def tearDownClass(cls):
        """Clean up test outputs but preserve reference data"""
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")

if __name__ == '__main__':
    unittest.main()
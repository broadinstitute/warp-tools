#!/usr/bin/env python3
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("mito_plot_filter.py")


class TestMitoPlotFilterScript(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SCRIPT_PATH.read_text(encoding="utf-8")

    def test_script_exists(self):
        self.assertTrue(SCRIPT_PATH.exists(), f"Missing script: {SCRIPT_PATH}")

    def test_uses_output_root_and_basename_for_outputs(self):
        expected_lines = [
            'output_prefix     = f"{output_root.rstrip(\'/\')}/{basename}"',
            'vcf_output_path             = f"{output_prefix}.filtered.vcf.bgz"',
            'sample_metadata_output_path = f"{output_prefix}.metadata.tsv"',
        ]
        for line in expected_lines:
            self.assertIn(line, self.source)

    def test_numt_warning_cast_to_int32_for_vcf(self):
        self.assertIn(
            "numt_fp_warning=hl.int32(mt_vcf.numt_fp_warning)",
            self.source,
            "VCF export must cast bool numt_fp_warning to int32",
        )

    def test_exports_vcf_with_tabix(self):
        self.assertIn("hl.export_vcf(", self.source)
        self.assertIn("vcf_output_path", self.source)
        self.assertIn("tabix=True", self.source)

    def test_cli_uses_output_root_and_basename(self):
        self.assertIn('parser.add_argument("--output-root", required=True', self.source)
        self.assertIn('parser.add_argument("--basename", required=True', self.source)

    def test_cli_and_entrypoint_present(self):
        self.assertIn("def parse_args():", self.source)
        self.assertIn("if __name__ == \"__main__\":", self.source)
        self.assertIn("main()", self.source)


if __name__ == "__main__":
    unittest.main(verbosity=2)

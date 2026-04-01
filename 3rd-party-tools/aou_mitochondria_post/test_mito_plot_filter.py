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

    def test_uses_output_base_for_local_outputs(self):
        expected_lines = [
            'vcf_local             = f"{output_base}.vcf.bgz"',
            'sample_metadata_local = f"{output_base}_metadata.tsv"',
            'variants_per_sample_svg            = f"{output_base}.variants_per_sample.svg"',
            'mito_cn_distribution_svg           = f"{output_base}.mito_cn_distribution.svg"',
            'variant_allele_frequency_svg       = f"{output_base}.variant_allele_frequency.svg"',
            'variant_af_and_allele_fraction_svg = f"{output_base}.variant_af_and_allele_fraction.svg"',
            'numt_fp_by_mtcn_svg                = f"{output_base}.numt_fp_by_mtcn.svg"',
            'haplogroup_heteroplasmy_svg        = f"{output_base}.haplogroup_heteroplasmy.svg"',
            'haplogroup_homoplasmy_svg          = f"{output_base}.haplogroup_homoplasmy.svg"',
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
        self.assertIn("hl.export_vcf(mt_vcf, vcf_local, tabix=True)", self.source)

    def test_cli_and_entrypoint_present(self):
        self.assertIn("def parse_args():", self.source)
        self.assertIn("if __name__ == \"__main__\":", self.source)
        self.assertIn("main()", self.source)


if __name__ == "__main__":
    unittest.main(verbosity=2)

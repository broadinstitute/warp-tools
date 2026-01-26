#!/usr/bin/env python3
"""
Minimal script to check for zero-valued metrics in a simple key/value or CSV file.

Usage:
  python3 test_zero_check.py path/to/metrics.csv

The script accepts files where each line is one of:
  metric_name: value
  metric_name\tvalue
  metric_name,value

It prints PASS if no problematic zero-valued metrics are found, otherwise
prints ERROR with the offending metric names and exits with code 1.
"""

import sys
import os


def _parse_numeric(value):
	if value is None:
		return None
	v = str(value).strip()
	if v == "":
		return None
	# allow percentages and thousands separators
	v = v.replace(',', '')
	if v.endswith('%'):
		v = v[:-1]
	try:
		return float(v)
	except ValueError:
		return None


def check_zero_values_in_library_metrics(file_path, file_name, skip_metrics=None):
	if skip_metrics is None:
		skip_metrics = [
			"reads_mapped_antisense_to_gene",
			"reads_mapped_confidently_to_intronic_regions",
			"percent_intronic_reads",
		]

	if not os.path.exists(file_path):
		print(f"ERROR: File not found: {file_path}")
		return False

	zero_metrics = []
	skipped_metrics = []

	try:
		with open(file_path, 'r') as f:
			for line in f:
				line = line.strip()
				if not line or line.startswith('#'):
					continue

				# Try simple key/value parsing (colon, tab, comma)
				if ':' in line and ',' not in line:
					key, value = line.split(':', 1)
				elif '\t' in line:
					parts = line.split('\t')
					if len(parts) >= 2:
						key, value = parts[0], parts[1]
					else:
						continue
				elif ',' in line:
					parts = line.split(',')
					if len(parts) >= 2:
						key, value = parts[0], parts[1]
					else:
						continue
				else:
					continue

				key = key.strip()
				value = value.strip()

				numeric_value = _parse_numeric(value)
				if numeric_value is None:
					continue

				if numeric_value == 0.0:
					if key in skip_metrics:
						skipped_metrics.append(key)
					else:
						zero_metrics.append(key)

		if skipped_metrics:
			print(f"INFO: {file_name} skipped zero-value metrics: {', '.join(skipped_metrics)}")

		if zero_metrics:
			print(f"ERROR: {file_name} has metrics with zero values: {', '.join(zero_metrics)}")
			return False
		else:
			print(f"PASS: {file_name} has no problematic metrics with zero values")
			return True

	except Exception as e:
		print(f"ERROR: Could not process {file_name}: {str(e)}")
		return False


if __name__ == "__main__":
	# Minimal CLI: take first arg as file path, else use the example name
	test_file = "nuclei_2k_mouse_example_1234_library_metrics.csv"
	if len(sys.argv) > 1:
		test_file = sys.argv[1]

	print(f"Testing library metrics file: {test_file}")
	print("=" * 60)

	result = check_zero_values_in_library_metrics(test_file, "Library Metrics")

	print("=" * 60)
	if result:
		print("✓ Test PASSED")
		sys.exit(0)
	else:
		print("✗ Test FAILED")
		sys.exit(1)


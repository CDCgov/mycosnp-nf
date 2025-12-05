#!/usr/bin/env python3

"""
QC report parser script to evaluate quality control metrics
against specified thresholds.
"""

import os
import sys
import argparse
import pandas as pd

def parse_thresholds(threshold_str):
    """Parse and validate the QC thresholds provided as a string."""
    try:
        parts = threshold_str.split(',')
        if len(parts) != 4:
            raise ValueError(
                "Thresholds must be 'GCrangePct:min-max,AvgQscore:#,RefLenCov:#,MeanCovDepth:#'"
            )

        # Extract the numbers for each threshold
        gc_range_str, coverage_str, depth_str, qscore_str = [part.split(':')[1] for part in parts]

        # Parse the GC range
        gc_range = list(map(float, gc_range_str.split('-')))
        if len(gc_range) != 2:
            raise ValueError("GC range must have exactly two values separated by a dash ('-').")

        # Parse the other thresholds
        coverage_threshold, depth_threshold, qscore_threshold = map(float, [coverage_str, depth_str, qscore_str])

        return gc_range, coverage_threshold, depth_threshold, qscore_threshold
    except ValueError as value_error:
        print(f"Error parsing thresholds: {value_error}")
        sys.exit(1)


def evaluate_qc(row, gc_range, coverage_threshold, depth_threshold, qscore_threshold):
    """Evaluate QC metrics for a single sample against the qc thresholds."""
    print(f"Evaluating QC for sample {row['Sample Name']}")
    if not gc_range[0] <= row['GC After Trimming'] <= gc_range[1]:
        print(f"GC After Trimming {row['GC After Trimming']} is out of range {gc_range}")
        return 'fail'
    if (
        row['Reference Length Coverage After Trimming'] >= coverage_threshold
        and row['Mean Coverage Depth'] >= depth_threshold
        and row['Average Q Score After Trimming'] >= qscore_threshold
    ):
        return 'PASS'

    print(f"Sample {row['Sample Name']} failed on thresholds.")
    return 'FAIL'


def main():
    """Main function to parse QC report, evaluate samples, and save results."""
    print("Starting the QC parser main script")

    parser = argparse.ArgumentParser(description='QC report parser')
    parser.add_argument(
        'file_path', help='Path to the QC report file'
    )
    parser.add_argument(
        '-qc_thresholds',
        required=True,
        help='QC thresholds format "GCrangePct:min-max,AvgQscore:#,RefLenCov:#,MeanCovDepth:#"',
    )
    args = parser.parse_args()

    # Check if the input file exists and is not empty
    if not os.path.isfile(args.file_path):
        print(f"Error: The file '{args.file_path}' does not exist.")
        sys.exit(1)

    if os.stat(args.file_path).st_size == 0:
        print(f"Error: The file '{args.file_path}' is empty.")
        sys.exit(1)

    # Parse and validate the thresholds
    gc_range, coverage_threshold, depth_threshold, qscore_threshold = parse_thresholds(
        args.qc_thresholds
    )

    try:
        data_frame = pd.read_csv(args.file_path, sep='\t')
        if data_frame.empty:
            raise ValueError("Input file is empty after parsing.")
        print("Successfully read the input file")
    except (ValueError, pd.errors.EmptyDataError) as read_error:
        print(f"Error reading the input file: {read_error}")
        sys.exit(1)

    try:
        data_frame = data_frame.replace('%', '', regex=True).apply(
            pd.to_numeric, errors='ignore'
        )
        print("Successfully converted columns to numeric")
    except (ValueError, TypeError) as conversion_error:
        print(f"Error during data conversion: {conversion_error}")
        sys.exit(1)

    try:
        data_frame['outcome'] = data_frame.apply(
            evaluate_qc,
            axis=1,
            gc_range=gc_range,
            coverage_threshold=coverage_threshold,
            depth_threshold=depth_threshold,
            qscore_threshold=qscore_threshold,
        )
        print("Successfully evaluated QC for all rows")
    except (ValueError, TypeError) as eval_error:
        print(f"Error during QC evaluation: {eval_error}")
        sys.exit(1)

    # Normalize column names: lowercase, replace spaces with underscores
    data_frame.columns = data_frame.columns.str.lower().str.replace(' ', '_')

    # Rename sample name column
    column_rename = {
        'sample_name': 'sample_id',
    }
    data_frame = data_frame.rename(columns=column_rename)

    # Print result
    try:
        data_frame.to_csv('combined_QC_results.csv', index=False)
        print("QC results saved to 'combined_QC_results.csv'.")
    except (IOError, OSError) as save_error:
        print(f"Error saving the results: {save_error}")
        sys.exit(1)


if __name__ == "__main__":
    main()

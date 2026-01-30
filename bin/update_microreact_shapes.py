#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
from datetime import datetime


def create_df(infile: str) -> pd.DataFrame:
    """return pandas df, columns in lowrecase"""
    return pd.read_csv(infile, low_memory=False, encoding='utf-8', encoding_errors='replace').rename(columns=str.lower)

#Checks if the sample names are in B####,B####,B#### format or B####, B####, B#### format
def format_test_samples(sample_list: list) -> list:
    """return list removing whitespace if found"""
    return [item.strip() for entry in sample_list for item in entry.split(',') if item.strip()]

def parse_date_columns(df: pd.DataFrame, date_column='sample_collection_date') -> pd.DataFrame:
    """parse date column and add day, month, year columns"""

    if date_column not in df.columns:
        raise KeyError(f"Column '{date_column}' not found in dataframe")

    def extract_date_parts(date_str):
        """extract day, month, year from date string, handle various formats"""
        if pd.isna(date_str) or date_str == '' or date_str == 'Not Provided':
            return pd.Series({'day': None, 'month': None, 'year': None})

        try:
            # Try MM/DD/YYYY format
            date_obj = datetime.strptime(str(date_str).strip(), '%m/%d/%Y')
            return pd.Series({
                'day': date_obj.day,
                'month': date_obj.month,
                'year': date_obj.year
            })
        except ValueError:
            try:
                # Try YYYY-MM-DD format
                date_obj = datetime.strptime(str(date_str).strip(), '%Y-%m-%d')
                return pd.Series({
                    'day': date_obj.day,
                    'month': date_obj.month,
                    'year': date_obj.year
                })
            except ValueError:
                # If parsing fails, return None
                return pd.Series({'day': None, 'month': None, 'year': None})

    # Apply date parsing
    date_parts = df[date_column].apply(extract_date_parts)
    df['day'] = date_parts['day'].astype('Int64')  # Use Int64 to handle NaN
    df['month'] = date_parts['month'].astype('Int64')
    df['year'] = date_parts['year'].astype('Int64')

    return df

def rename_colour_column(df: pd.DataFrame) -> pd.DataFrame:
    """rename clade_color to clade_colour if it exists"""
    if 'clade_color' in df.columns:
        df = df.rename(columns={'clade_color': 'clade_colour'})
    return df

def validate_sample_ids(samplesheet_file: str, metadata_df: pd.DataFrame) -> None:
    """
    Validate that sample IDs in samplesheet match those in metadata.
    Exits with error code 1 if validation fails.

    Args:
        samplesheet_file: Path to the input samplesheet CSV
        metadata_df: DataFrame containing metadata with sample_id column
    """

    # Read samplesheet
    try:
        df_samplesheet = create_df(samplesheet_file)
    except Exception as e:
        print(f"\nERROR: Failed to read samplesheet file '{samplesheet_file}'", file=sys.stderr)
        print(f"Details: {str(e)}\n", file=sys.stderr)
        sys.exit(1)

    # Check for required columns
    if 'sample_id' not in metadata_df.columns:
        print("\nERROR: Metadata file must contain a 'sample_id' column", file=sys.stderr)
        print(f"Available columns: {', '.join(metadata_df.columns)}\n", file=sys.stderr)
        sys.exit(1)

    # Extract sample IDs (remove empty/null values)
    samplesheet_samples = set(df_samplesheet['sample'].dropna().astype(str).str.strip())
    metadata_samples = set(metadata_df['sample_id'].dropna().astype(str).str.strip())

    # Remove empty strings
    samplesheet_samples = {s for s in samplesheet_samples if s}
    metadata_samples = {s for s in metadata_samples if s}

    # Check if samplesheet samples are in metadata
    missing_in_metadata = samplesheet_samples - metadata_samples
    extra_in_metadata = metadata_samples - samplesheet_samples

    # Report validation results only if there are mismatches
    if missing_in_metadata or extra_in_metadata:
        print("\nERROR: Sample ID validation failed", file=sys.stderr)

        if missing_in_metadata:
            print(f"The following {len(missing_in_metadata)} sample(s) in samplesheet are NOT found in metadata:\n", file=sys.stderr)
            for sample in sorted(missing_in_metadata):
                print(f"  {sample}", file=sys.stderr)
            print("", file=sys.stderr)

        if extra_in_metadata:
            print(f"The following {len(extra_in_metadata)} sample(s) in metadata are NOT found in samplesheet:\n", file=sys.stderr)
            for sample in sorted(extra_in_metadata):
                print(f"  {sample}", file=sys.stderr)
            print("", file=sys.stderr)

        print("Sample IDs in samplesheet and metadata must match exactly.\n", file=sys.stderr)
        sys.exit(1)

def add_geolocation(df_meta: pd.DataFrame, df_geo: pd.DataFrame, shared_column='state_country') -> pd.DataFrame:
    """check df for correct state_country code, convert if needed, and merge left df and df_geo"""

    if shared_column not in df_meta.columns:
        raise KeyError(f"Column '{shared_column}' not found in df")
    if shared_column not in df_geo.columns:
        raise KeyError(f"Column '{shared_column}' not found in df_geo")

    missing_cells = set(df_meta[shared_column]) - set(df_geo[shared_column])

    """
    logic:
    check which values in df_meta[shared_column] are not in df_geo[shared_column]
        missing value is scanned over the entire df_geo data to find which row matches value
            find which row has the value
                if the value exist and is only a single row
                    rewrite the value correctly using the df_geo[shared_column]
                else
                    do nothing and that values remains the same and data will not be merged from df_geo for that row in df_meta
    """
    if missing_cells:
        # convert string to lowercase in case user
        df_geo_lower = df_geo.astype(str).map(str.lower)
        for i in missing_cells:
            i_lower = i.lower()
            row = df_geo_lower[df_geo_lower.eq(i_lower).any(axis=1)]
            if not row.empty and len(row[shared_column]) == 1: # len > 1 suggests geolocation file has mutiple rows with identical cell values
                valid_cell = row[shared_column].iloc[0].upper()
                df_meta.loc[df_meta[shared_column] == i, shared_column] = valid_cell

    # merge both dfs on shared_column
    merged_df = pd.merge(df_meta, df_geo, on=shared_column, how='left', suffixes=('', '_geo'))
    merged_df = merged_df.drop(columns=merged_df.filter(regex="_geo$").columns)

    columns = {
        'state_latitude':'latitude__state_country',
        'state_longitude':'longitude__state_country',
        'region_lat':'latitude__arln_region',
        'region_long':'longitude__arln_region'
    }

    return merged_df.rename(columns=columns)

def add_shapes(df: pd.DataFrame, sample_list: list) -> pd.DataFrame:
    """vectorize to add values to given columns"""
    df.insert(7,'state_country__Shape','Circle')
    # isolate samples
    if sample_list:
        mask = df['sample_id'].isin(sample_list)
        df.loc[mask, 'state_country__Shape'] = 'Star'
        df.loc[mask, 'state_country'] += '_new'
    return df

#Creates the new metadata microreact csv
def create_new_csv(df: pd.DataFrame, output='microreact_metadata.csv') -> None:
    df.to_csv(output, index=False)

def opts() -> argparse.Namespace:
    """"adding dedicated parse funtion for better readability"""
    description = """
        Script intakes two csv files, one containing metadata for microreact_metadata (positional) and one containing geolocaitons (-g)
          the sample names or -t flag applies a new column with shapes and appends "_new" to the values in state_country
    """

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        dest='microreact_metadata',
        type=str,
        help='CSV file with microreact metadata'
    )
    parser.add_argument(
        '-t','--test_samples',
        nargs='*',
        type=str,
        help='Comma separated list of test samples',
        dest="test_samples",
        metavar="sample_id"
    )
    parser.add_argument(
        '-g','--geolocation',
        type=str,
        help='CSV file with state_country centroid lat and long',
        dest='geolocation',
        metavar="file"
    )
    parser.add_argument(
        '-s','--samplesheet',
        type=str,
        help='Input samplesheet CSV file to validate sample IDs against metadata',
        dest='samplesheet',
        metavar="file",
        default=None
    )

    return parser.parse_args()

def main():

    args = opts()
    # params
    microreact_metadata = args.microreact_metadata
    test_samples = format_test_samples(sample_list=args.test_samples) if args.test_samples else args.test_samples
    geolocation = args.geolocation

    # import files
    df_metadata_raw = create_df(infile=microreact_metadata)
    df_geolocation_raw = create_df(infile=geolocation)

    # validate sample IDs if samplesheet provided
    if args.samplesheet:
        validate_sample_ids(samplesheet_file=args.samplesheet, metadata_df=df_metadata_raw)

    # rename color column and parse dates
    df_metadata_raw = rename_colour_column(df=df_metadata_raw)
    df_metadata_raw = parse_date_columns(df=df_metadata_raw)

    # merge df's and add shapes column
    metadata_shape_geo_df = add_geolocation(df_meta=df_metadata_raw, df_geo=df_geolocation_raw)
    metadata_shape_geo_df = add_shapes(df=metadata_shape_geo_df, sample_list=test_samples)

    # write output
    create_new_csv(df=metadata_shape_geo_df)

if __name__ == "__main__":

    main()

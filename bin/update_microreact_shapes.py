#!/usr/bin/env python3

import argparse
import pandas as pd


def create_df(infile: str) -> pd.DataFrame:
    """return pandas df, columns in lowrecase"""
    return pd.read_csv(infile, low_memory=False, encoding='utf-8', encoding_errors='replace').rename(columns=str.lower)

#Checks if the sample names are in B####,B####,B#### format or B####, B####, B#### format
def format_test_samples(sample_list: list) -> list:
    """return list removing whitespace if found"""
    return [item.strip() for entry in sample_list for item in entry.split(',') if item.strip()]

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

    # merge df's and add shapes column
    metadata_shape_geo_df = add_geolocation(df_meta=df_metadata_raw, df_geo=df_geolocation_raw)
    metadata_shape_geo_df = add_shapes(df=metadata_shape_geo_df, sample_list=test_samples)

    # write output
    create_new_csv(df=metadata_shape_geo_df)

if __name__ == "__main__":

    main()

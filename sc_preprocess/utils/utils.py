import os
import sys
import argparse

import pandas as pd

from sc_preprocess.utils.custom_logger import custom_logger


def validate_cores(value):
    """Validate that cores is either a positive integer or 'all'."""
    if value.lower() == 'all':
        return value
    try:
        cores = int(value)
        if cores < 1:
            raise argparse.ArgumentTypeError("cores must be a positive integer or 'all'")
        return str(cores)
    except ValueError:
        raise argparse.ArgumentTypeError("cores must be a positive integer or 'all'")


# Function to apply suffix if specified
def get_directories(config):
    """Get directory paths with optional suffix applied"""
    dirs = {}
    base_dirs = config.get("directories")
    suffix = config.get("directories_suffix", "none")
    
    for key, value in base_dirs.items():
        if suffix != "none" and suffix != "":
            # Apply suffix to directory path
            dirs[key] = f"{value}_{suffix}"
        else:
            # Use directory path as-is
            dirs[key] = value
    
    return dirs

def has_underscore(s: str) -> bool:
    """
    Return True if the string contains at least one underscore, else False.
    """
    return "_" in s


def sanity_check_libraries_list_tsv(filepath, log_file=None, expected_columns=None, path_column="fastqs", file_extension=None, validate_csv_contents=True):
    """
    Check the format of the libraries list TSV file.
    Args:
        filepath (str): Path to the TSV file.
        log_file (str, optional): Path to the log file. Defaults to None.
        expected_columns (set, optional): Expected column names. Defaults to {"batch", "capture", "sample", "fastqs"}.
        path_column (str, optional): Name of the column containing file paths. Defaults to "fastqs".
        file_extension (str, optional): Required file extension (e.g., ".csv"). If None, checks for directory existence.
    Returns:
        pd.DataFrame: Valid dataframe or exits on error
    """
    # Set default expected columns if not provided
    if expected_columns is None:
        expected_columns = {"batch", "capture", "sample", "fastqs"}

    # Read the file
    try:
        df = pd.read_csv(filepath, sep="\t")
        df.columns = df.columns.str.strip()

    except Exception as e:
        custom_logger.error(f"Pandas could not read your libraries filepath here: '{filepath}'. "
                            f"This was the error: {e}")
        sys.exit(1)

    # Validate headers
    actual_columns = set(df.columns)
    if expected_columns != actual_columns:
        custom_logger.error(f"Expected columns {expected_columns}, found {actual_columns}")
        sys.exit(1)

    # Check for duplicate batch+capture combinations
    dup_cols = [c for c in ["batch", "capture"] if c in df.columns]
    if len(dup_cols) == 2:
        duplicates = df[df.duplicated(subset=dup_cols, keep=False)]
        if not duplicates.empty:
            for _, dup_row in duplicates.iterrows():
                custom_logger.error(f"Duplicate batch+capture combination found: batch='{dup_row['batch']}', capture='{dup_row['capture']}'. "
                                    f"Each batch+capture pair must appear only once.")
            sys.exit(1)

    valid = True
    for idx, row in df.iterrows():
        if has_underscore(str(row["batch"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Batch name '{row['batch']}' cannot contain underscores.")
            valid = False

        if "capture" in expected_columns and has_underscore(str(row["capture"])):
            custom_logger.error(f"Row {idx + 1} in '{filepath}': Capture name '{row['capture']}' cannot contain underscores.")
            valid = False
            
        file_path = row[path_column]
        error_location = f"Row {idx + 2} in '{filepath}'"
        if not isinstance(file_path, str):
            custom_logger.error(f"{error_location}: Path is not a string: {file_path}")
            valid = False
            continue

        # Handle different path formats (fastqs vs CSV)
        if path_column == "fastqs":
            # For GEX/ATAC: Handle comma-separated fastq directories
            DELIMITERS = [",", ";", ":", " "]
            # First, check for disallowed delimiters
            for delim in DELIMITERS[1:]:  # skip comma (allowed)
                if delim in file_path:
                    custom_logger.error(f"{error_location}: Path contains invalid delimiter '{delim}'. Only comma-separated paths are allowed.")
                    sys.exit(1)

            # Split the paths by comma
            suspected_paths = [path.strip() for path in file_path.split(",") if path.strip()]
            # Convert relative paths to absolute and validate existence
            absolute_paths = []
            for path in suspected_paths:
                if not os.path.isabs(path):
                    # Convert relative path to absolute
                    abs_path = os.path.abspath(path)
                    custom_logger.warning(f"{error_location}: Converting relative path '{path}' to absolute path '{abs_path}'")
                    path = abs_path
                
                if not os.path.exists(path):
                    custom_logger.error(f"{error_location}: Path does not exist: {path}")
                    valid = False
                else:
                    absolute_paths.append(path)
            
            # Update the dataframe with absolute paths
            if absolute_paths:
                df.at[idx, path_column] = ",".join(absolute_paths)
        else:
            # For ARC CSV or other single file formats
            # Check file extension if specified
            if file_extension and not file_path.endswith(file_extension):
                custom_logger.error(f"{error_location}: The path '{file_path}' does not end with '{file_extension}'. Please ensure all paths have the proper file extension.")
                valid = False
                continue
            
            # Convert relative path to absolute
            if not os.path.isabs(file_path):
                abs_path = os.path.abspath(file_path)
                custom_logger.warning(f"{error_location}: Converting relative path '{file_path}' to absolute path '{abs_path}'")
                file_path = abs_path
                df.at[idx, path_column] = abs_path
            
            if not os.path.exists(file_path):
                custom_logger.error(f"{error_location}: File does not exist: {file_path}")
                valid = False
            elif file_extension == ".csv" and validate_csv_contents:
                # For ARC per-capture CSVs: validate fastqs paths inside the CSV
                try:
                    lib_csv = pd.read_csv(file_path)
                    if "fastqs" not in lib_csv.columns:
                        custom_logger.error(f"{error_location}: '{file_path}' is missing a 'fastqs' column.")
                        valid = False
                    else:
                        csv_modified = False
                        expanded_rows = []
                        for csv_idx, csv_row in lib_csv.iterrows():
                            csv_loc = f"Row {csv_idx + 2} in '{file_path}'"
                            raw = str(csv_row["fastqs"]).strip()
                            paths = [p.strip() for p in raw.split(",") if p.strip()]
                            if len(paths) > 1:
                                csv_modified = True
                            for fastqs_path in paths:
                                if not os.path.isabs(fastqs_path):
                                    abs_fastqs = os.path.abspath(fastqs_path)
                                    custom_logger.warning(f"{csv_loc}: Converting relative path '{fastqs_path}' to absolute path '{abs_fastqs}'")
                                    fastqs_path = abs_fastqs
                                    csv_modified = True
                                if not os.path.exists(fastqs_path):
                                    custom_logger.error(f"{csv_loc}: Path does not exist: {fastqs_path}")
                                    valid = False
                                else:
                                    new_row = csv_row.copy()
                                    new_row["fastqs"] = fastqs_path
                                    expanded_rows.append(new_row)
                        if csv_modified and expanded_rows:
                            pd.DataFrame(expanded_rows).to_csv(file_path, index=False)
                except Exception as e:
                    custom_logger.error(f"{error_location}: Could not read '{file_path}': {e}")
                    valid = False

    if valid:
        custom_logger.info(f"{filepath} file format is valid.")
        return df
    else:
        custom_logger.error(f"Some errors were found in {filepath}.")
        sys.exit(1)
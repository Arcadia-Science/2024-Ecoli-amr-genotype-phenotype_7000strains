import sys

import pandas as pd


def extract_mpileup_info(input_file, output_file):
    # Dynamically determine the number of initial comment lines to skip
    with open(input_file) as file:
        skip_count = sum(1 for line in file if line.startswith("#"))

    # Read the mpileup file into a pandas DataFrame
    # Assuming the mpileup file has no header and columns are separated by tabs
    col_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    mpileup_data = pd.read_csv(
        input_file, sep="\t", header=None, names=col_names, skiprows=skip_count
    )

    # Extract the DP value from the Info column
    mpileup_data["DP"] = mpileup_data["INFO"].str.extract(r"DP=(\d+)")

    # Convert the DP column to numeric,
    # setting errors to 'coerce' which replaces non-convertible values with NaN
    mpileup_data["DP"] = pd.to_numeric(mpileup_data["DP"], errors="coerce")

    # Replace NaN values with -1 (or another placeholder if -1 is not suitable)
    mpileup_data["DP"].fillna(-1, inplace=True)

    # Convert to integer
    mpileup_data["DP"] = mpileup_data["DP"].astype(int)

    # Select only the Chromosome, Position, and DP columns
    extracted_data = mpileup_data[["CHROM", "POS", "DP"]]

    # Save the extracted data to a tab-delimited file
    extracted_data.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file> <output_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_mpileup_info(input_file, output_file)

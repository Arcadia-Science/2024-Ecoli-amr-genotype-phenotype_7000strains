import logging

import numpy as np
import pandas as pd

# Set up logging to keep track of the process
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Step 1: Read Locus list from File1
## To be able to obtain DP for each ECOR strain at each position in the pangenome,
## we provide the script with a txt file that contains all the positions of pangenomes
logging.info("Reading loci from List of LOCI file")

with open("dataset_generation/data/dp_threshold/index_loci_pangenome_good.txt") as f:
    loci = [line.strip() for line in f]
loci_index = {locus: i for i, locus in enumerate(loci)}  # Mapping loci to indices

# Step 2: Initialize NA1 array
##Initialize an array to store the DP information.
##Each row is a position/nucleotide in the pangenome
##The posistion is identified by the locus name and position within the locus)
##Each column represents one strain of the ECOR collection
##We provide a .txt file that defines which text files containing the DP information need to be used

logging.info("Reading file paths from the list of samples")
with open("dataset_generation/data/dp_threshold/list_ecor_txtfiles.txt") as f:
    file_list = [line.strip() for line in f]
num_loci = len(loci)
num_files = len(file_list)
na1 = np.zeros((num_loci, num_files), dtype=int)

# Step 3: Process each data file
for file_idx, file_path in enumerate(file_list):
    logging.info(f"Processing file {file_idx + 1}/{num_files}: {file_path}")
    data = pd.read_csv(file_path, sep=r"\s+", usecols=["CHROM", "POS", "DP"])
    data["LOCUS"] = data["CHROM"].astype(str) + "_" + data["POS"].astype(str)

    # Update NA1 array
    for _, row in data.iterrows():
        locus = row["LOCUS"]
        if locus in loci_index:
            na1[loci_index[locus], file_idx] = row["DP"]

# Include the LOCUS information
enhanced_na1 = np.column_stack((loci, na1))

# Save the final array to a text file
output_file = "ecor72_array.txt"
header = "\t".join(["LOCUS"] + file_list)
logging.info(f"Saving the final Ecor72 array to {output_file}")
np.savetxt(output_file, enhanced_na1, fmt="%s", delimiter="\t", header=header, comments="")

logging.info("Process completed successfully.")

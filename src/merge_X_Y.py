import csv
import sys
from pathlib import Path
from paths import PATHS

# merge the two csv files into one based on the first column (pdb_id)
def merge_csv(y_labels, x_labels):
    with open(y_labels, "r") as y_file, open(x_labels, "r") as x_file:
        y_reader = csv.reader(y_file)
        x_reader = csv.reader(x_file)

        # Read headers and detect key columns.
        x_header = next(x_reader)
        y_header = next(y_reader)
        x_key_col = "pdb_id" if "pdb_id" in x_header else ("ID" if "ID" in x_header else x_header[0])
        y_key_col = "pdb_id" if "pdb_id" in y_header else ("ID" if "ID" in y_header else y_header[0])
        x_key_idx = x_header.index(x_key_col)
        y_key_idx = y_header.index(y_key_col)

        # Create a new header for the merged CSV
        merged_header = y_header + [col for col in x_header if col != x_key_col]

        # Create a dictionary to store the merged data
        merged_data = {}

        # Read the y_labels and store them in the merged_data dictionary
        for row in y_reader:
            if not row:
                continue
            pdb_id = row[y_key_idx]
            merged_data[pdb_id] = row

        # Read the x_labels and merge them with the corresponding y_labels
        for row in x_reader:
            if not row:
                continue
            pdb_id = row[x_key_idx]
            if pdb_id in merged_data:
                merged_data[pdb_id] += [val for idx, val in enumerate(row) if idx != x_key_idx]

    return merged_header, list(merged_data.values())

if __name__ == "__main__":
    y_labels = Path(sys.argv[1]) if len(sys.argv) > 1 else PATHS.inp_pdb_dir.joinpath("ylabels.csv")
    x_labels = PATHS.features_dir.joinpath("custom_features_merged_all.csv")

    merged_header, merged_data = merge_csv(y_labels, x_labels)

    # Write the merged data to a new CSV file
    merged_csv_path = PATHS.features_dir.joinpath("merged_features.csv")
    with open(merged_csv_path, "w", newline="") as merged_file:
        writer = csv.writer(merged_file)
        writer.writerow(merged_header)
        writer.writerows(merged_data)

    print(f"Merged CSV file created at: {merged_csv_path}")


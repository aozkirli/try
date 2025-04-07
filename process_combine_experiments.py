import os
import glob
import re
import shutil
from collections import Counter

import pandas as pd

# --- Main Script ---

# Get the full path of the current script.
# In interactive environments __file__ might not be defined so we fallback to the current working directory.
try:
    script_path = os.path.abspath(__file__)
except NameError:
    script_path = os.path.abspath(os.getcwd())

script_dir = os.path.dirname(script_path)
print("Script directory:", script_dir)

# Define the data directory relative to the script directory.
data_dir = os.path.abspath(os.path.join(script_dir, 'data'))
print("Data directory (assumed):", data_dir)

# Remove existing 'studies' folder if it exists, then create a new one.
studies_path = os.path.join(data_dir, 'studies')
if os.path.exists(studies_path):
    shutil.rmtree(studies_path)
os.makedirs(studies_path, exist_ok=True)

# Change directory to the 'experiments' folder inside the data directory.
experiments_dir = os.path.join(data_dir, 'experiments')
os.chdir(experiments_dir)
print("Experiments Directory:", os.getcwd())

# List all CSV files in the 'experiments' folder.
csv_files = glob.glob("*.csv")
csv_files = [str(f) for f in csv_files]

# Extract experiment names from file names using a regex.
pattern = re.compile(r'^(.*?_20\d{2}\w?)(?=[_.])')
expnames = [pattern.search(f).group(1) if pattern.search(f) else "" for f in csv_files]

# Tabulate occurrences of each unique experiment name.
exp_counter = Counter(expnames)
summary = list(exp_counter.items())

print("Experiment Summary:")
for name, count in summary:
    print(f"  {name}: {count}")

# Loop over each unique experiment name.
unique_expnames = [name for name, count in summary]
for k, target_name in enumerate(unique_expnames, start=1):
    os.chdir(experiments_dir)
    df_accum = pd.DataFrame()

    # Find indices of all files that match the current experiment name.
    indices = [i for i, ename in enumerate(expnames) if ename == target_name]

    # Process each matching file.
    for i, idx in enumerate(indices, start=1):
        filename = csv_files[idx]
        print(f"Processing file: {filename}")
        df = pd.read_csv(filename,delimiter=';')
        df['expnum'] = i  # Add experiment number column.
        df_accum = pd.concat([df_accum, df], ignore_index=False)

    # Use the last processed file name for renaming.
    rename = csv_files[indices[-1]]
    # Force the .csv extension if it's missing:
    if not rename.lower().endswith('.csv'):
        rename += '.csv'
    
    # Remove '_ExpX' tags and spaces from the file name.
    count_target = exp_counter[target_name]
    for jj in range(1, count_target + 1):
        rename = rename.replace(f"_Exp{jj}", "").replace(" ", "")
    
    # Add an index number to the file name with zero-padding.
    num = f"{k:02d}"
    rename = f"{num}_{rename}"
    
    # Save the cumulative DataFrame as a CSV file with ';' as the delimiter.
    csv_save_path = os.path.join(studies_path, rename)
    df_accum.to_csv(csv_save_path, sep=';', index=False)
    print("Saved:", csv_save_path)  

# Return to the original script directory.
os.chdir(script_dir)
print("Returned to script directory:", os.getcwd())

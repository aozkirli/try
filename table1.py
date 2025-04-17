import pandas as pd
import numpy as np
import os

# Read the master table CSV
tbl = pd.read_csv('data/datasets/SD_ma_master_table.csv', delimiter=';')
# Ensure numeric columns are floats
num_cols = tbl.select_dtypes(include='number').columns
tbl[num_cols] = tbl[num_cols].astype(np.float64)

# Initialize list to collect summary rows
rows = []
variable = 'error_ori_deb'

# Unique dataset codes
datasets = tbl['studynum'].unique()

# Loop over each dataset
for i in datasets:
    tmp = tbl[tbl['studynum'] == i]

    # Compute fields
    ID = str(int(i))      
    Study = tmp['study'].iloc[0]
    Stimulus = tmp['stimulus'].iloc[0]
    Datasets = tmp['codenum'].nunique()
    N = tmp['obsid'].nunique()
    AllTrials = len(tmp)
    CleanTrials = tmp.loc[tmp[variable].notna() & (tmp['delta'].abs() <= 90)].shape[0]

    rows.append({
        'ID': ID,
        'Study': Study,
        'Stimulus': Stimulus,
        'Datasets': Datasets,
        'N': N,
        'AllTrials': AllTrials,
        'CleanTrials': CleanTrials
    })

# Create summary DataFrame
report_table = pd.DataFrame(rows)

# Add a total/summary row using concat
total_row = {
    'ID': '',
    'Study': '',
    'Stimulus': 'Total',
    'Datasets': report_table['Datasets'].sum(),
    'N': report_table['N'].sum(),
    'AllTrials': report_table['AllTrials'].sum(),
    'CleanTrials': report_table['CleanTrials'].sum()
}
report_table = pd.concat([report_table, pd.DataFrame([total_row])], ignore_index=True)

# Ensure output directory exists
out_dir = 'tables'
os.makedirs(out_dir, exist_ok=True)

# Save to CSV
out_path = os.path.join(out_dir, 'summary_studies.csv')
report_table.to_csv(out_path, index=False)

print(f"Saved summary table to {out_path}")

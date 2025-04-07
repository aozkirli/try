import os
import glob
import numpy as np
import pandas as pd

# =============================================================================
# Helper functions
# =============================================================================

def nanunique(series):
    """Return the unique values in a pandas Series, ignoring NaNs."""
    return series.dropna().unique()

def is_outlier_quartiles(series):
    """Return a boolean Series marking outliers using the 1.5*IQR rule."""
    q1 = series.quantile(0.25)
    q3 = series.quantile(0.75)
    iqr = q3 - q1
    return (series < q1 - 1.5*iqr) | (series > q3 + 1.5*iqr)

def zscore(series):
    """Return a z-scored version of the pandas Series."""
    return (series - series.mean()) / series.std()

# =============================================================================
# Main Script
# =============================================================================

# Determine script directory (if __file__ is not defined—for example, in a notebook—use cwd)
try:
    script_path = os.path.abspath(__file__)
except NameError:
    script_path = os.path.abspath(os.getcwd())
script_dir = os.path.dirname(script_path)
print("Script directory:", script_dir)

# Change directory to the 'studies' folder
studies_path = os.path.abspath(os.path.join(script_dir, 'data', 'studies'))
os.chdir(studies_path)
print("Current studies directory:", os.getcwd())

# List all .mat files in the directory
csvfiles = glob.glob("*.csv")
csvfiles = sorted(csvfiles)  # sort if needed
print(csvfiles)

# For this translation, we assume each .mat file corresponds to one study.
# (In MATLAB, study_list = unique(cellstr(csvfiles)); so here we simply use the file names.)
study_list = np.unique(csvfiles)

# Standardized variable names to keep
variables = ['obs', 'theta', 'resp', 'delta', 'error', 'cond', 
             'rt', 'block', 'study', 'experiment', 'stimulus', 'expnum']

print("\nGenerating Master table\n")
master_tbl = pd.DataFrame()  # initialize empty master table

# Loop over each study MAT file to load and process its data
for i, csvfile in enumerate(study_list, start=1):
    print(f"Processing study file: {csvfile}")
    try:
        tbl = pd.read_table(csvfile, delimiter=';')
    except Exception as e:
        print(f"Error loading {csvfile}: {e}")
        continue

    # Add missing columns if needed:
    if 'cond' not in tbl.columns:
        tbl['cond'] = 1
    if 'block' not in tbl.columns:
        tbl['block'] = 1
    if 'rt' not in tbl.columns:
        tbl['rt'] = 9  # default value

    # Keep only standardized variables that exist in the DataFrame
    cols_to_keep = [v for v in variables if v in tbl.columns]
    tbl = tbl[cols_to_keep].copy()

    # Add 'stimtype' (flag for 'Orientation') and 'studynum'
    # (Assuming tbl['stimulus'] is string; convert to lowercase for comparison)
    print("DataFrame columns:", tbl.columns.tolist())
    tbl['stimtype'] = tbl['stimulus'].astype(str).str.lower() == 'orientation'
    tbl['studynum'] = i

    # Append to master table
    master_tbl = pd.concat([master_tbl, tbl], ignore_index=True)

    print(f"\nMaster table ready with {i} studies combined.\n")

# =============================================================================
# Recoding observer IDs and handling dataset specifics
# =============================================================================

print("\nRecoding each dataset and incremental observer ID...\n")
master = master_tbl.copy()
id_study = master['study'].unique()
n = len(id_study)
# Preallocate an array for number of subjects per study (assume up to 5 experiments per study)
n_subjects = np.full((n, 5), np.nan)
master_coded = pd.DataFrame()
count = 0

for i, study in enumerate(id_study):
    tbl = master[master['study'] == study].copy()
    id_expe = tbl['experiment'].unique()
    n_exp = len(id_expe)
    
    for j, expe in enumerate(id_expe, start=1):
        tbl_e = tbl[tbl['experiment'] == expe].copy()
        # Count unique observers (ignoring NaN)
        unique_obs = pd.Series(tbl_e['obs']).dropna().unique()
        n_subjects[i, j-1] = len(unique_obs)
        
        # Fix subject numbering if necessary:
        if len(unique_obs) > 0 and len(unique_obs) != tbl_e['obs'].max():
            mapping = {old: new for new, old in enumerate(np.sort(unique_obs), start=1)}
            tbl_e['obs'] = tbl_e['obs'].map(mapping)
        
        id_cond = nanunique(tbl_e['cond'])
        n_cond = len(id_cond)
        
        for k, cond in enumerate(id_cond, start=1):
            tbl_e_k = tbl_e[tbl_e['cond'] == cond].copy()
            count += 1
            
            # Generate a new code for this dataset
            numcode = f"{count:02d}"
            name = f"{numcode} {study}"
            if n_exp > 1:
                name += f" E{j}"
            if n_cond > 1:
                name += f" C{k}"
            tbl_e_k['code'] = name  # assign code for all rows
            
            # Update 'obsid' with incremental IDs
            if master_coded.empty:
                tbl_e_k['obsid'] = tbl_e_k['obs']
            else:
                # Compare current study and experiment with last row of master_coded
                last_row = master_coded.iloc[-1]
                if (tbl_e_k['study'].iloc[0] == last_row['study'] and 
                    tbl_e_k['expnum'].iloc[0] == last_row['expnum']):
                    new_obsid = master_coded['obsid'].max() - n_subjects[i, j-1] + tbl_e_k['obs']
                    tbl_e_k['obsid'] = new_obsid
                else:
                    print(f"Dataset: {tbl_e_k['study'].iloc[0]}, {tbl_e_k['experiment'].iloc[0]}, "
                          f"n_subjects: {n_subjects[i, j-1]}, "
                          f"prev max obsid: {master_coded['obsid'].max()}")
                    tbl_e_k['obsid'] = master_coded['obsid'].max() + tbl_e_k['obs']
            
            master_coded = pd.concat([master_coded, tbl_e_k], ignore_index=True)

master = master_coded.copy()
print("\nRecoding done.\n")

# =============================================================================
# Error Preprocessing
# =============================================================================

nbasis = 6  # Number of basis functions for bias removal
studies_unique = master['study'].unique()
new_tbl = pd.DataFrame()

for s in studies_unique:
    tbl_study = master[master['study'] == s].copy()
    # Round theta and delta to integers
    tbl_study['theta'] = tbl_study['theta'].round()
    tbl_study['delta'] = tbl_study['delta'].round()
    
    for expe in tbl_study['experiment'].unique():
        tbl_experiment = tbl_study[tbl_study['experiment'] == expe].copy()
        
        for cond in tbl_experiment['cond'].unique():
            tbl_condition = tbl_experiment[tbl_experiment['cond'] == cond].copy()
            
            for subj in tbl_condition['obs'].unique():
                tbl_subject = tbl_condition[tbl_condition['obs'] == subj].copy()
                
                # Remove improbable data: errors > 90°
                out = tbl_subject['error'].abs() > 90
                erroriqr = tbl_subject['error'].copy()
                errorsd = tbl_subject['error'].copy()
                erroriqr[out] = np.nan
                errorsd[out] = np.nan
                
                # Remove outliers based on quartiles (IQR method)
                out_q = is_outlier_quartiles(erroriqr)
                erroriqr[out_q] = np.nan
                
                # Remove outliers based on 3 standard deviations
                mean_abs = np.nanmean(errorsd.abs())
                std_val = np.nanstd(errorsd)
                threshold = mean_abs + 3 * std_val
                out_sd = errorsd.abs() > threshold
                errorsd[out_sd] = np.nan
                
                # Remove outliers based on reaction time (> 10s)
                out_rt = tbl_subject['rt'] > 10
                erroriqr[out_rt] = np.nan
                errorsd[out_rt] = np.nan
                
                tbl_subject['outliers_q'] = out_q
                tbl_subject['outliers_sd'] = out_sd
                tbl_subject['outliers_rt'] = out_rt
                
                tbl_subject['error_iqr'] = erroriqr
                tbl_subject['errorsd'] = errorsd

                # Orientation Bias Correction
                if str(tbl_subject['stimulus'].iloc[0]).lower() == 'orientation':
                    theta360 = tbl_subject['theta'] * 2
                    tbl_subject['theta_cent'] = tbl_subject['theta'] - 90
                else:
                    theta360 = tbl_subject['theta']
                    tbl_subject['theta_cent'] = tbl_subject['theta'] - 180

                N = len(tbl_subject)
                Xsin = np.empty((N, nbasis))
                Xcos = np.empty((N, nbasis))
                for b in range(1, nbasis+1):
                    Xsin[:, b-1] = np.sin(np.deg2rad(theta360 * b))
                    Xcos[:, b-1] = np.cos(np.deg2rad(theta360 * b))
                X = np.hstack([np.ones((N, 1)), Xsin, Xcos])
                # Identify rows with NaN in erroriqr or in the first basis column (after intercept)
                nanout = np.isnan(tbl_subject['error_iqr']) | np.isnan(X[:, 1])
                # Fit regression using pseudo-inverse on non-NaN rows
                beta = np.linalg.pinv(X[~nanout, :]) @ tbl_subject['error_iqr'][~nanout]
                predicted = X @ beta
                peak = np.nanmax(predicted)
                tbl_subject['stim_bias_peak'] = peak  # replicate peak value for all rows
                tbl_subject['stim_bias_peak'] = np.full((N,), peak)
                tbl_subject['error_ori_deb'] = tbl_subject['error_iqr'] - predicted

                # Serial Dependence Bias Correction
                if str(tbl_subject['stimulus'].iloc[0]).lower() == 'orientation':
                    delta360 = tbl_subject['delta'] * 2
                else:
                    delta360 = tbl_subject['delta']
                Xsin_delta = np.empty((N, nbasis))
                Xcos_delta = np.empty((N, nbasis))
                for b in range(1, nbasis+1):
                    Xsin_delta[:, b-1] = np.sin(np.deg2rad(delta360 * b))
                    Xcos_delta[:, b-1] = np.cos(np.deg2rad(delta360 * b))
                X_delta = np.hstack([Xsin_delta, Xcos_delta])
                error_ori_deb = tbl_subject['error_ori_deb']
                nanout_delta = np.isnan(error_ori_deb) | np.isnan(X_delta[:, 0])
                beta_delta = np.linalg.pinv(X_delta[~nanout_delta, :]) @ error_ori_deb[~nanout_delta]
                predicted_delta = X_delta @ beta_delta
                tbl_subject['error_ori_deb_sd_deb'] = tbl_subject['error_ori_deb'] - predicted_delta
                peak_delta = np.nanmax(predicted_delta)
                tbl_subject['sd_bias_peak'] = np.full((N,), peak_delta)

                # Normalized versions of error measures (z-score normalization)
                tbl_subject['error_norm'] = zscore(tbl_subject['error'])
                tbl_subject['error_iqr_norm'] = zscore(tbl_subject['error_iqr'])
                tbl_subject['error_ori_deb_norm'] = zscore(tbl_subject['error_ori_deb'])
                tbl_subject['error_ori_deb_sd_deb_norm'] = zscore(tbl_subject['error_ori_deb_sd_deb'])

                # Create bins based on delta values
                tbl_subject['bin'] = np.nan
                tbl_subject.loc[tbl_subject['delta'].abs() <= 10, 'bin'] = 1
                tbl_subject.loc[(tbl_subject['delta'].abs() >= 40) & (tbl_subject['delta'].abs() <= 50), 'bin'] = 2
                tbl_subject.loc[(tbl_subject['delta'].abs() >= 80) & (tbl_subject['delta'].abs() <= 90), 'bin'] = 3

                # Append preprocessed subject data to new_tbl
                new_tbl = pd.concat([new_tbl, tbl_subject], ignore_index=True)

# Add a numeric version of the code (convert to categorical codes)
new_tbl['codenum'] = new_tbl['code'].astype('category').cat.codes.astype(float)
master = new_tbl.copy()

# =============================================================================
# Save final datasets
# =============================================================================

# Change directory to the datasets folder
datasets_path = os.path.abspath(os.path.join(script_dir, 'data', 'datasets'))
os.chdir(datasets_path)
print("Saving final datasets to:", os.getcwd())

# Save as CSV (using semicolon as delimiter)
master.to_csv('SD_ma_master_table.csv', sep=';', index=False)

# Make adjustments: set theta==180 to 0 and delta==-90 to 90
tbl = master.copy()
tbl.loc[tbl['theta'] == 180, 'theta'] = 0
tbl.loc[tbl['delta'] == -90, 'delta'] = 90

# Return to the original script directory
os.chdir(script_dir)
print("Returned to script directory:", os.getcwd())

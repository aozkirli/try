import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load preprocessed data 
tb = pd.read_table('data/datasets/SD_ma_master_table.csv', delimiter=';')
num_cols = tb.select_dtypes(include='number').columns
tb[num_cols] = tb[num_cols].astype(np.float64)

# Initialize parameters and containers
error_var = 'error_ori_deb'
bin_size = 11
degrees = [2, 3]

tbl_scatter_list = []
best_fit_all = []
best_polynomial_degree = []
mvav_scatter_aggregated = []
mvav_bias_aggregated = []

# Loop over study numbers
for studynum in tb['studynum'].unique():
    tbl_i = tb[tb['studynum'] == studynum].copy()
    tbl_i = tbl_i[np.abs(tbl_i['delta']) <= 90]

    # If sd_deb not already removed, flip sign for negative delta
    if 'sd_deb' not in error_var:
        mask = (tbl_i['delta'] < 0) & (tbl_i['delta'] != -90)
        tbl_i.loc[mask, error_var] *= -1

    # Loop over experiments
    for expnum in tbl_i['expnum'].unique():
        tbl_i_k = tbl_i[tbl_i['expnum'] == expnum]

        # Loop over conditions
        for cond in tbl_i_k['cond'].unique():
            tbl_i_k_j = tbl_i_k[tbl_i_k['cond'] == cond]

            # Loop over observers
            for obs in tbl_i_k_j['obs'].unique():
                tbl_i_k_j_o = tbl_i_k_j[tbl_i_k_j['obs'] == obs]

                # Compute grouped stats
                abs_delta = np.abs(tbl_i_k_j_o['delta'])
                scatt = tbl_i_k_j_o.groupby(abs_delta)[error_var].std()
                biass = tbl_i_k_j_o.groupby(abs_delta)[error_var].mean()
                counts = tbl_i_k_j_o.groupby(abs_delta).size()

                # Prepare fit arrays
                delta_vals = scatt.index.values
                mask_nz = (scatt.fillna(0) != 0)
                delta_fit = delta_vals[mask_nz]
                scatt_fit = scatt[mask_nz].values
                bias_fit = biass[mask_nz].values
                w = counts[mask_nz].values
                weights = w / w.max()

                # Skip if too few points
                if len(delta_fit) < min(degrees) + 1:
                    continue

                # BIC-driven model selection
                BIC_values = []
                for deg in degrees:
                    p = np.polyfit(delta_fit, scatt_fit, deg, w=np.sqrt(weights))
                    y_pred = np.polyval(p, delta_fit)
                    resid = scatt_fit - y_pred
                    SSR = np.sum((np.sqrt(weights) * resid) ** 2) / weights.sum()
                    num_params = deg + 1
                    n = len(delta_fit)
                    BIC = n * np.log(SSR / n) + num_params * np.log(n)
                    BIC_values.append(BIC)

               # Sort BIC and select best degree within +2 of minimum
                sorted_indices = np.argsort(BIC_values)
                min_bic = BIC_values[sorted_indices[0]]
                eligible = [idx for idx in sorted_indices if BIC_values[idx] < min_bic + 2]
                best_deg = degrees[min(eligible)]
                best_polynomial_degree.append(best_deg)

                # Final fit & prediction
                p_best = np.polyfit(delta_fit, scatt_fit, best_deg, w=np.sqrt(weights))
                fit_pred = np.polyval(p_best, np.arange(0, 91))
                best_fit_all.append(fit_pred)

                # Discard if non-positive or NaN at key indices
                key_idx = [0, 45, 90]
                if np.any(fit_pred[key_idx] <= 0) or np.any(np.isnan(fit_pred[key_idx])):
                    continue

                # Helper to compute 1D weighted moving average
                def weighted_mvav(values, wts, window):
                    ser_num = pd.Series(values * wts)
                    ser_den = pd.Series(wts)
                    mvnum = ser_num.rolling(window, center=True, min_periods=1).sum()
                    mvden = ser_den.rolling(window, center=True, min_periods=1).sum()
                    return (mvnum / mvden).values

                # Build full-symmetric series (extend by reflection)
                def reflect_and_aggregate(vals, cnts, sign):
                    base = np.full(91, np.nan)
                    base[np.array(delta_fit, int)] = vals
                    refl = sign*base[-2:0:-1]
                    ext_vals = np.tile(np.concatenate([refl, base]), 3)

                    weights = np.full(91, np.nan)
                    weights[np.array(delta_fit, int)] = cnts
                    refl = weights[-2:0:-1]
                    ext_cnts = np.tile(np.concatenate([refl, weights]), 3)
                    return ext_vals, ext_cnts

                # Bias moving average
                mv_vals, mv_wts = reflect_and_aggregate(bias_fit, w,-1)
                mv_bias = weighted_mvav(mv_vals, mv_wts, bin_size)
                mvav_bias_aggregated.append(mv_bias[269 + np.arange(91)])

                # Scatter moving average
                mv_vals, mv_wts = reflect_and_aggregate(scatt_fit, w,1)
                mv_scatt = weighted_mvav(mv_vals, mv_wts, bin_size)
                mvav_scatter_aggregated.append(mv_scatt[269 + np.arange(91)])

                # Compute bin-wise std of error for 'bin' categories, skipping NaN bins
                bin_std = tbl_i_k_j_o.dropna(subset=['bin']).groupby('bin')[error_var].std()
            
                # Store scatter table rows
                rows = pd.DataFrame({
                    'obsid': tbl_i_k_j_o['obsid'].iloc[0],
                    'codenum': tbl_i_k_j_o['codenum'].iloc[0],
                    'SI': ['iso', 'mid', 'ortho'],
                    'ES': fit_pred[key_idx],
                    'bin_scatter': bin_std,
                    'ntrials': len(w)
                })
                tbl_scatter_list.append(rows)


# Consolidate and save
tbl_scatter = pd.concat(tbl_scatter_list, ignore_index=True)
best_fit_array = np.column_stack(best_fit_all)
mvav_bias_arr = np.column_stack(mvav_bias_aggregated)
mvav_scatt_arr = np.column_stack(mvav_scatter_aggregated)

# Save to disk
tbl_scatter.to_csv('tbl_scatter.csv', sep=';', float_format='%.16f', index=False)
np.savetxt('best_fit_all.csv', best_fit_array, delimiter=';')
np.savetxt('mvav_bias_aggregated.csv', mvav_bias_arr, delimiter=';')
np.savetxt('mvav_scatter_aggregated.csv', mvav_scatt_arr, delimiter=';')
print('Results saved: tbl_scatter.csv, best_fit_all.csv, mvav_bias_aggregated.csv, mvav_scatter_aggregated.csv')


import numpy as np
import pandas as pd
from scipy.stats import t, chi2
import matplotlib.pyplot as plt
import os

# Load MATLAB results
mvav_bias = pd.read_csv('mvav_bias_aggregated.csv',delimiter=';',dtype=float, header=None)
mvav_scatter = pd.read_csv('mvav_scatter_aggregated.csv',delimiter=';',dtype=float, header=None)
best_fit = pd.read_csv('best_fit_all.csv',delimiter=';',dtype=float, header=None)

# Create figure with two side-by-side panels
fig, axes = plt.subplots(1, 2, figsize=(20, 10), sharey=False)

# Panel 1: bias
plt.sca(axes[0])
x_plot = np.arange(0, 91)
y_plot = mvav_bias.mean(axis=1).values
mv_n_sm = np.isfinite(mvav_bias).sum(axis=1)
mv_std_sm = mvav_bias.std(axis=1).values
dfree = mv_n_sm - 1
alpha = 0.05
tcrit = t.ppf(1 - alpha / 2, dfree)
ci_lower = y_plot - tcrit * mv_std_sm / np.sqrt(mv_n_sm)
ci_upper = y_plot + tcrit * mv_std_sm / np.sqrt(mv_n_sm)

p, = plt.plot(x_plot, y_plot, color='k', linewidth=1.5)
p.set_linewidth(4)
p.set_linestyle(':')
f = plt.fill_between(x_plot, ci_upper, ci_lower, color='k', alpha=0.2)
axes[0].set_xlabel('|Δ| (°)', fontsize=20)
axes[0].set_ylabel('Bias (°)', fontsize=20)
axes[0].set_title('Data', fontsize=20)
axes[0].set_xticks([0, 45, 90])
axes[0].set_xticklabels(['iso', 'mid', 'ortho'])
axes[0].grid(True)
axes[0].tick_params(labelsize=20, width=1.9)
plt.ylim(-.5,2)

# Panel 2: scatter + fit
plt.sca(axes[1])
y_plot = mvav_scatter.mean(axis=1).values
mv_n_sm = np.isfinite(mvav_scatter).sum(axis=1)
mv_std_sm = mvav_scatter.std(axis=1).values
dfree = mv_n_sm - 1
tcrit = t.ppf(1 - alpha / 2, dfree)
ci_lower = y_plot - tcrit * mv_std_sm / np.sqrt(mv_n_sm)
ci_upper = y_plot + tcrit * mv_std_sm / np.sqrt(mv_n_sm)

p, = plt.plot(x_plot, y_plot, color='k', linewidth=1.5)
p.set_linewidth(4)
p.set_linestyle(':')
f = plt.fill_between(x_plot, ci_upper, ci_lower, color='k', alpha=0.2)

# poly fit mean and ci
y_plot = best_fit.mean(axis=1).values
mv_n_sm = np.isfinite(best_fit).sum(axis=1)
mv_std_sm = best_fit.std(axis=1).values
dfree = mv_n_sm - 1
tcrit = t.ppf(1 - alpha / 2, dfree)
ci_lower = y_plot - tcrit * mv_std_sm / np.sqrt(mv_n_sm)
ci_upper = y_plot + tcrit * mv_std_sm / np.sqrt(mv_n_sm)
p, = plt.plot(x_plot, y_plot, color='c', linewidth=1.5)
p.set_linewidth(4)
f = plt.fill_between(x_plot, ci_upper, ci_lower, color='c', alpha=0.2)

axes[1].set_xlabel('|Δ| (°)', fontsize=20)
axes[1].set_ylabel('Error Scatter (°)', fontsize=20)
axes[1].set_title('Data vs. Polynomial Fit', fontsize=20)
axes[1].set_xticks([0, 45, 90])
axes[1].set_xticklabels(['iso', 'mid', 'ortho'], fontsize=20)
axes[1].grid(True)
plt.ylim(8.8,10.6)
axes[1].tick_params(labelsize=20, width=1.9)

# Save figure
os.makedirs('figures', exist_ok=True)
out_pdf = os.path.join('figures', 'main_results_polyfit.pdf')
fig.savefig(out_pdf, format='pdf', transparent=True)
# Also save a raster version
out_png = os.path.join('figures', 'main_results_polyfit.png')
fig.savefig(out_png, dpi=300, transparent=True)
plt.close(fig)

print(f"Saved figure to {out_pdf} and {out_png}")
"""Fig. 6: monthly burned cropland and CO2, whole basin, 2020.
Run analyse_comparison_gfed.py first (writes results_basin_monthly.csv)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
m = pd.read_csv('results_basin_monthly.csv', index_col=0)
lab = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
area = {'Sentinel-2 + VIIRS (this study)': m.s2_ha, 'MODIS MCD64A1': m.modis_ha, 'FireCCI51': m.firecci_ha}
co2 = {'Sentinel-2 + VIIRS (this study)': m.s2_tCO2, 'MODIS MCD64A1': m.modis_tCO2, 'FireCCI51': m.firecci_tCO2,
       'GFED4.1s agriculture': m.gfed_tCO2}
col = {'Sentinel-2 + VIIRS (this study)': '#d95f0e', 'MODIS MCD64A1': '#3b6e9c', 'FireCCI51': '#7fb1d6',
       'GFED4.1s agriculture': '#7a7a7a'}
plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(2, 1, figsize=(12, 8.4), dpi=250)
for a, d, yl, t in [(ax[0], area, 'Burned cropland (10$^3$ ha)', '(a) Burned cropland area'),
                    (ax[1], co2, 'CO$_2$ (kt)', '(b) CO$_2$ emissions')]:
    n = 4; w = 0.8 / n; x = np.arange(12)
    for i, (k, v) in enumerate(d.items()):
        a.bar(x - 0.4 + w / 2 + i * w, np.asarray(v) / 1e3, w, label=k, color=col[k], edgecolor='white', lw=0.5)
    a.set_xticks(x); a.set_xticklabels(lab); a.set_ylabel(yl); a.set_title(t, loc='left', fontweight='bold')
    a.spines[['top', 'right']].set_visible(False); a.grid(axis='y', color='#ddd', lw=0.6); a.set_axisbelow(True)
h, l = ax[1].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig('fig6_comparison.png', bbox_inches='tight')

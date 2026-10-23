"""Tier 2 comparison, whole basin, January-December 2020.
Sums the 12 monthly CSVs exported by 02_basin_S2_vs_MODIS_FireCCI.js by crop and month,
computes overlaps and the MCD64A1 + Sentinel-2-only combination, and extracts GFED4.1s
agricultural CO2 for the same area and months.
Usage: python analyse_comparison_gfed.py results_basin_S2_MODIS_FireCCI_m*.csv GFED4.1s_2020_beta.hdf5
Writes results_basin_comparison_summary.json and results_basin_monthly.csv."""
import sys, json, numpy as np, pandas as pd
NAMES = {1: 'Corn', 2: 'Cotton', 3: 'Rice', 4: 'Sorghum', 5: 'Soybeans', 10: 'Peanuts', 13: 'Popcorn',
         24: 'Winter wheat', 28: 'Oats', 29: 'Millet', 45: 'Sugarcane', 46: 'Sweet potatoes'}
csvs = [a for a in sys.argv[1:] if a.lower().endswith('.csv')]
h5 = [a for a in sys.argv[1:] if a.lower().endswith(('.hdf5', '.h5'))]
df = pd.concat([pd.read_csv(f) for f in csvs])
df = df[df.crop.isin(NAMES)].copy(); df['name'] = df.crop.map(NAMES)
# t CO2 per ha of each crop, to convert overlap areas
df['f'] = (df.s2_tCO2 / df.s2_ha).where(df.s2_ha > 0, df.modis_tCO2 / df.modis_ha)
df['f'] = df.f.where(df.f.notna(), df.firecci_tCO2 / df.firecci_ha).fillna(0)
df['s2_only_ha'] = df.s2_ha - df.s2_and_modis_ha
df['s2_only_tCO2'] = df.s2_tCO2 - df.s2_and_modis_ha * df.f
df['modis_plus_s2only_ha'] = df.modis_ha + df.s2_only_ha
df['modis_plus_s2only_tCO2'] = df.modis_tCO2 + df.s2_only_tCO2
num = [c for c in df.columns if c.endswith('_ha') or c.endswith('_tCO2')]
bym = df.groupby('month')[num].sum(); byc = df.groupby('name')[num].sum().sort_values('s2_tCO2', ascending=False)
tot = df[num].sum(); tot['cropland_ha'] = df[df.month == df.month.min()].cropland_ha.sum()
pd.set_option('display.width', 250); pd.set_option('display.float_format', '{:,.0f}'.format)
print(bym); print(byc); print(tot)
out = {'by_month': bym.round(1).to_dict('index'), 'by_crop': byc.round(1).to_dict('index'), 'total': tot.round(1).to_dict()}
out['s2_found_by_modis'] = tot.s2_and_modis_ha / tot.s2_ha
out['modis_found_by_s2'] = tot.s2_and_modis_ha / tot.modis_ha
if h5:
    import h5py
    from matplotlib.path import Path
    aoi = np.array(json.load(open('aoi_coords.json'))[0])
    f = h5py.File(h5[0], 'r')
    lat = f['lat'][:]; lon = f['lon'][:]; area = f['ancill/grid_cell_area'][:]
    # fraction of each 0.25 deg cell inside the AOI, via 10 x 10 sub-sampling
    ii = np.where((lat[:, 0] > 29) & (lat[:, 0] < 39))[0]; jj = np.where((lon[0] > -95) & (lon[0] < -86))[0]
    p = Path(aoi); frac = np.zeros(lat.shape); off = (np.arange(10) + 0.5) / 10 * 0.25 - 0.125
    for i in ii:
        for j in jj:
            pts = np.c_[(lon[i, j] + 0 * off[:, None] + off[None, :]).ravel(), (lat[i, j] + off[:, None] + 0 * off[None, :]).ravel()]
            frac[i, j] = p.contains_points(pts).mean()
    EF_CO2_AGRI = 1585.0  # g CO2 / kg DM, GFED4 agricultural (Akagi et al. 2011)
    g = {}
    for m in range(1, 13):
        grp = f['emissions/%02d' % m]
        dm_agri = grp['DM'][:] * grp['partitioning/DM_AGRI'][:] * area * frac / 1e3   # t DM
        g[m] = {'DM_agri_t': float(dm_agri.sum()), 'CO2_agri_t': float(dm_agri.sum() * EF_CO2_AGRI / 1e3)}
    out['gfed'] = g; out['gfed_cells'] = int((frac > 0).sum())
    bym['gfed_tCO2'] = [g[m]['CO2_agri_t'] for m in bym.index]
    print(json.dumps(g, indent=1))
bym.to_csv('results_basin_monthly.csv')
json.dump(out, open('results_basin_comparison_summary.json', 'w'), indent=1, default=float)

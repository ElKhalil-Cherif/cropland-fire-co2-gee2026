import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
fig, ax = plt.subplots(figsize=(11, 12.5), dpi=200); ax.set_xlim(0, 110); ax.set_ylim(0, 125); ax.axis('off')
C = {'in': '#dff3f7', 't1': '#fde7c8', 't2': '#e6e1f5', 'em': '#f9d0e6', 'out': '#d8f0d2'}
def box(x, y, w, h, t, c, bold=False):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle='round,pad=0.4,rounding_size=1.2', fc=C[c], ec='#555', lw=1))
    ax.text(x, y, t, ha='center', va='center', fontsize=9.5, fontweight='bold' if bold else 'normal', wrap=True)
def arr(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='-|>', color='#333', lw=1.1))
def group(x, y, w, h, t):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1.5', fc='none', ec='#999', ls='--', lw=1))
    ax.text(x + 1.5, y + h - 1.8, t, fontsize=10, fontweight='bold', va='top')
# inputs
group(2, 104, 106, 19, 'Data inputs (Google Earth Engine)')
box(20, 111, 30, 7, 'USDA CDL 2020\n(30 m crop types, 12 crops)', 'in')
box(55, 111, 30, 7, 'FAOSTAT 2020 crop yields\n(United States, kg ha$^{-1}$)', 'in')
box(90, 111, 30, 7, 'Residue parameters (Table 1)\nRPR, DM, CF, EF$_{CO_2}$', 'in')
# tier 1
group(2, 47, 44, 55, 'Tier 1: whole basin, Jan-Dec 2020')
box(24, 93, 38, 7, 'MODIS MCD64A1 burn date\n(500 m, monthly)', 't1')
box(24, 80, 38, 7, 'Cropland mask (CDL)\n120 m grid, EPSG:5070', 't1')
box(24, 67, 38, 7, 'Monthly burn masks and\nannual union (each pixel once)', 't1')
box(24, 54, 38, 7, 'VIIRS 375 m consistency check\n(point-in-pixel and 375 m)', 't1')
for a, b in [(93, 80), (80, 67), (67, 54)]: arr(24, a - 3.6, 24, b + 3.6)
# tier 2
group(50, 31, 58, 71, 'Tier 2: Testing Zone (a), Sep-Nov 2020')
box(79, 93, 52, 7, 'Sentinel-2 L2A (<20% cloud), SCL cloud/\nshadow/water mask with 60 m buffer', 't2')
box(79, 82, 52, 7, 'Pre-fire median (30 days) and post-fire\nlowest-NBR composite (month + 15 days)', 't2')
box(79, 71, 52, 7, 'Change indices: dNBR, dMIRBI, dBAIS2,\ndNDTI, dSAVI', 't2')
box(79, 60, 52, 7, 'Analysis mask: VIIRS fires buffered 3 km\nand CDL cropland', 't2')
box(79, 49, 52, 7, 'Otsu threshold per index and S2 tile;\nburned if at least 4 of 5 indices agree', 't2')
box(79, 38, 52, 7, 'Remove patches < 3000 m$^2$\n30 m monthly burn masks', 't2')
for a, b in [(93, 82), (82, 71), (71, 60), (60, 49), (49, 38)]: arr(79, a - 3.6, 79, b + 3.6)
# emission
box(55, 22, 70, 8, 'Crop emission coefficient (t CO$_2$ ha$^{-1}$) = Yield x RPR x DM x CF x EF / 10$^6$\nE = burned area x coefficient, per pixel and crop', 'em', True)
arr(48, 103.6, 48, 26.5)
arr(24, 50.4, 36, 26.5); arr(79, 34.4, 70, 26.5)
# outputs
box(19, 7, 34, 8, 'Annual and monthly CO$_2$\nper crop, emission map', 'out')
box(55, 7, 32, 8, 'Monte Carlo uncertainty\n(10,000 draws)', 'out')
box(91, 7, 34, 8, 'Intercomparison in Zone (a):\nSentinel-2 vs MCD64A1 vs FireCCI51', 'out')
arr(45, 17.8, 22, 11.2); arr(55, 17.8, 55, 11.2); arr(65, 17.8, 88, 11.2)
fig.savefig('fig2_flowchart.png', bbox_inches='tight')

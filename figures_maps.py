import json, math, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from PIL import Image
NE = 'naturalearth/'   # Natural Earth GeoJSON files (see README)
AOI = [[-89.833008,37.308108],[-91.606201,35.715921],[-91.861084,35.019116],[-92.243408,34.664959],[-92.243408,34.266222],[-91.518311,33.655409],[-91.584229,33.215125],[-92.221436,32.430002],[-91.836914,31.834296],[-92.485108,31.084225],[-92.287354,30.792170],[-93.667236,30.109015],[-91.965910,29.738164],[-90.477905,29.661294],[-90.322998,29.920179],[-90.427368,30.108103],[-90.744324,30.152648],[-90.793213,30.705001],[-91.363349,30.721486],[-91.520508,31.040631],[-90.853638,32.144007],[-89.925293,33.279201],[-89.914307,34.229610],[-89.474854,35.104721],[-88.486084,35.599556],[-88.272949,36.469881],[-88.290528,37.049043],[-88.244385,37.358198],[-88.002686,37.749851],[-89.004639,37.913901],[-89.917603,37.904376],[-90.038452,37.541083],[-89.833008,37.308108]]
EXT = [-95.2, -86.3, 28.6, 38.6]
def rings(geom):
    t, c = geom['type'], geom['coordinates']
    if t == 'Polygon': return c
    if t == 'MultiPolygon': return [r for p in c for r in p]
    if t == 'LineString': return [c]
    if t == 'MultiLineString': return c
    return []
def load(name): return json.load(open(NE + name + '.geojson'))['features']
def base(ax, ext=EXT, land=True):
    ax.set_facecolor('#d6e9f5' if land else 'white')
    if land:
        for f in load('ne_50m_land'):
            for r in rings(f['geometry']):
                xs, ys = zip(*r); ax.fill(xs, ys, color='#f4f1ea', lw=0, zorder=0)
    for f in load('ne_50m_admin_1_states_provinces_lakes'):
        if f['properties'].get('admin') != 'United States of America': continue
        for r in rings(f['geometry']):
            xs, ys = zip(*r); ax.plot(xs, ys, color='#9a9a9a', lw=0.6, zorder=1)
    for f in load('ne_10m_lakes'):
        for r in rings(f['geometry']):
            xs, ys = zip(*r)
            if ext[0] - 2 < min(xs) < ext[1] + 2 and ext[2] - 2 < min(ys) < ext[3] + 2:
                ax.fill(xs, ys, color='#d6e9f5', lw=0, zorder=1)
    for f in load('ne_10m_rivers_lake_centerlines'):
        n = f['properties'].get('name')
        if n in ('Mississippi', 'Arkansas', 'Red', 'Ohio', 'Yazoo', 'Atchafalaya', 'White', 'Ouachita'):
            for r in rings(f['geometry']):
                xs, ys = zip(*r); ax.plot(xs, ys, color='#3b83c4', lw=1.6 if n == 'Mississippi' else 0.8, zorder=2)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_aspect(1 / math.cos(math.radians((ext[2] + ext[3]) / 2)))
    ax.set_xlabel('Longitude (°)'); ax.set_ylabel('Latitude (°)')
def aoi(ax, lw=1.6):
    xs, ys = zip(*AOI); ax.plot(xs, ys, color='#111', lw=lw, zorder=5)
def north(ax, x, y):
    ax.annotate('N', xy=(x, y), xytext=(x, y - 0.75), ha='center', va='top', fontsize=13, fontweight='bold',
                arrowprops=dict(arrowstyle='-|>', lw=1.8, color='#111'), zorder=8)
def scalebar(ax, x, y, km=100):
    lat = y; dlon = km / (111.32 * math.cos(math.radians(lat)))
    ax.add_patch(Rectangle((x, y), dlon / 2, 0.08, fc='#111', ec='#111', zorder=8))
    ax.add_patch(Rectangle((x + dlon / 2, y), dlon / 2, 0.08, fc='white', ec='#111', zorder=8))
    for v, t in [(0, '0'), (dlon / 2, str(km // 2)), (dlon, str(km))]:
        ax.text(x + v, y + 0.15, t, ha='center', va='bottom', fontsize=9, zorder=8)
    ax.text(x + dlon + 0.12, y + 0.04, 'km', ha='left', va='center', fontsize=9, zorder=8)
CITIES = {'Memphis': (-90.05, 35.15), 'Little Rock': (-92.29, 34.75), 'Jackson': (-90.18, 32.30), 'Baton Rouge': (-91.15, 30.45),
          'Lafayette': (-92.02, 30.22), 'New Orleans': (-90.07, 29.95), 'St. Louis': (-90.20, 38.63), 'Greenville': (-91.06, 33.41),
          'Monroe': (-92.12, 32.51), 'Shreveport': (-93.75, 32.52)}
def cities(ax, names):
    for n in names:
        x, y = CITIES[n]; ax.plot(x, y, 'o', ms=4, mfc='white', mec='#111', zorder=9)
        ax.text(x + 0.08, y + 0.06, n, fontsize=8.5, zorder=9, color='#222')

def fig1(out):
    plt.rcParams.update({'font.size': 11})
    fig, ax = plt.subplots(figsize=(8.6, 10), dpi=220)
    base(ax); aoi(ax)
    cities(ax, ['Memphis', 'Little Rock', 'Jackson', 'Baton Rouge', 'Lafayette', 'New Orleans', 'Greenville', 'Monroe', 'Shreveport'])
    for n, (x, y) in {'ARKANSAS': (-93.3, 35.0), 'MISSISSIPPI': (-89.6, 32.9), 'LOUISIANA': (-92.9, 31.2), 'TENNESSEE': (-87.9, 35.9),
                      'MISSOURI': (-91.9, 36.9), 'ILLINOIS': (-89.0, 38.35), 'KENTUCKY': (-87.4, 37.15), 'ALABAMA': (-87.2, 32.6)}.items():
        ax.text(x, y, n, fontsize=8, color='#777', ha='center', style='italic', zorder=3)
    ax.text(-88.2, 29.0, 'Gulf of Mexico', fontsize=10, color='#3b6e9c', style='italic', ha='center')
    north(ax, -86.9, 38.4); scalebar(ax, -94.9, 28.85)
    h = [Line2D([0], [0], color='#111', lw=1.6, label='Study area boundary'),
         Line2D([0], [0], color='#3b83c4', lw=1.6, label='Mississippi River and tributaries'),
         Line2D([0], [0], color='#9a9a9a', lw=0.6, label='State boundaries')]
    ax.legend(handles=h, loc='center right', fontsize=8.5, framealpha=0.95, bbox_to_anchor=(1.0, 0.24))
    ins = ax.inset_axes([-95.15, 36.75, 2.5, 1.8], transform=ax.transData)
    for f in load('ne_50m_admin_0_countries'):
        if f['properties'].get('ADMIN') in ('United States of America', 'Canada', 'Mexico'):
            for r in rings(f['geometry']):
                xs, ys = zip(*r); ins.fill(xs, ys, color='#e6e1d6', ec='#888', lw=0.3)
    xs, ys = zip(*AOI); ins.fill(xs, ys, color='#d62728', lw=0)
    ins.set_xlim(-125, -66); ins.set_ylim(24, 50); ins.set_aspect('auto'); ins.set_xticks([]); ins.set_yticks([])
    ins.set_facecolor('#d6e9f5')
    fig.savefig(out, bbox_inches='tight')

def fig3(thumb, out):
    plt.rcParams.update({'font.size': 12})
    pal = ['#fff7bc', '#fec44f', '#fe9929', '#d95f0e', '#993404']
    xs, ys = zip(*AOI); ext = [min(xs), max(xs), min(ys), max(ys)]
    fig, ax = plt.subplots(figsize=(8.6, 11.5), dpi=220)
    e = [-94.6, -87.4, 29.2, 38.4]
    base(ax, e); 
    img = Image.open(thumb).convert('RGBA')
    ax.imshow(img, extent=ext, origin='upper', interpolation='nearest', zorder=4)
    aoi(ax, 1.2)
    ax.set_xlim(e[0], e[1]); ax.set_ylim(e[2], e[3])
    north(ax, -94.1, 38.2); scalebar(ax, -94.4, 29.45)
    ax.set_title('Annual cropland-fire CO$_2$ emissions, 2020')
    h = [Line2D([0], [0], color='#111', lw=1.2, label='Study area boundary'),
         Line2D([0], [0], color='#3b83c4', lw=1.6, label='Mississippi River')]
    ax.legend(handles=h, loc='center right', fontsize=9.5, framealpha=0.95, bbox_to_anchor=(1.0, 0.36))
    sm = plt.cm.ScalarMappable(cmap=LinearSegmentedColormap.from_list('e', pal), norm=Normalize(0, 2000))
    cb = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.03)
    cb.set_label('t CO$_2$ per km$^2$ of burned cropland (> 2000 shown at maximum)')
    fig.savefig(out, bbox_inches='tight')

def fig3pair(thumb, cropthumb, out):
    plt.rcParams.update({'font.size': 12})
    pal = ['#fff7bc', '#fec44f', '#fe9929', '#d95f0e', '#993404']
    xs, ys = zip(*AOI); ext = [min(xs), max(xs), min(ys), max(ys)]
    e = [-94.6, -87.4, 29.2, 38.4]
    fig, axs = plt.subplots(1, 2, figsize=(16, 11.2), dpi=220, gridspec_kw={'wspace': 0.08})
    for ax, img, title in [(axs[0], thumb, '(a) CO$_2$ emission intensity'), (axs[1], cropthumb, '(b) Crop type of burned cropland')]:
        base(ax, e)
        ax.imshow(Image.open(img).convert('RGBA'), extent=ext, origin='upper', interpolation='nearest', zorder=4)
        aoi(ax, 1.2)
        ax.set_xlim(e[0], e[1]); ax.set_ylim(e[2], e[3])
        north(ax, -94.1, 38.2); scalebar(ax, -94.4, 29.45)
        ax.set_title(title, loc='left', fontsize=14, fontweight='bold')
    axs[1].set_ylabel('')
    ref = [Line2D([0], [0], color='#111', lw=1.2, label='Study area boundary'),
           Line2D([0], [0], color='#3b83c4', lw=1.6, label='Mississippi River')]
    axs[0].legend(handles=ref, loc='center right', fontsize=10, framealpha=0.95, bbox_to_anchor=(1.0, 0.36))
    crops = [('Soybeans', '#8c2d04'), ('Rice', '#ff9f1c'), ('Sugarcane', '#6a3d9a'), ('Corn', '#ffd23f'),
             ('Cotton', '#1f78b4'), ('Other crops', '#9e9e9e')]
    axs[1].legend(handles=[Patch(fc=c, ec='none', label=n) for n, c in crops] + ref[:1], loc='center right',
                  fontsize=10, framealpha=0.95, bbox_to_anchor=(1.0, 0.36), title='Burned crop (2020)')
    sm = plt.cm.ScalarMappable(cmap=LinearSegmentedColormap.from_list('e', pal), norm=Normalize(0, 2000))
    cax = axs[0].inset_axes([0.90, 0.06, 0.03, 0.18])
    cb = fig.colorbar(sm, cax=cax); cb.ax.tick_params(labelsize=9, left=True, labelleft=True, right=False, labelright=False)
    cb.ax.yaxis.set_label_position('left')
    cb.set_label('t CO$_2$ km$^{-2}$ of burned cropland\n(> 2000 shown at maximum)', fontsize=9)
    fig.savefig(out, bbox_inches='tight')

if __name__ == '__main__':
    # map_thumbnail.png: MAP_THUMBNAIL_URL of script 01; crop_map.png: CROP_MAP_URL of script 03
    fig1('fig1_study_area.png'); fig3pair('map_thumbnail.png', 'crop_map.png', 'fig3_map.png')

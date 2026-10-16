"""Same-area comparison, NORTH_HALF (AOI north of 33.75 N), Sep-Nov 2020.
Usage: python3 analyse_comparison_gfed.py CF2020_v2_NORTH_HALF_m09.csv ..._m10.csv ..._m11.csv [GFED4.1s_2020_beta.hdf5]
GFED4.1s files: https://www.globalfiredata.org
Writes summary_north.json and prints tables."""
import sys, json, numpy as np, pandas as pd
NAMES={1:'Corn',2:'Cotton',3:'Rice',4:'Sorghum',5:'Soybeans',10:'Peanuts',13:'Popcorn',
       24:'Winter wheat',28:'Oats',29:'Millet',45:'Sugarcane',46:'Sweet potatoes'}
csvs=[a for a in sys.argv[1:] if a.lower().endswith('.csv')]
h5=[a for a in sys.argv[1:] if a.lower().endswith(('.hdf5','.h5'))]
df=pd.concat([pd.read_csv(f) for f in csvs])
df['crop']=df['crop'].astype(int); df['name']=df.crop.map(NAMES)
num=[c for c in df.columns if c.endswith('_ha') or c.endswith('_tCO2')]
out={}
bym=df.groupby('month')[num].sum(); byc=df.groupby('name')[num].sum().sort_values('s2_tCO2',ascending=False)
tot=df[num].sum()
pd.set_option('display.width',250); pd.set_option('display.float_format','{:,.0f}'.format)
print(bym); print(byc); print(tot)
out['by_month']=bym.round(1).to_dict('index'); out['by_crop']=byc.round(1).to_dict('index'); out['total']=tot.round(1).to_dict()
r=lambda a,b: a/b if b else float('nan')
out['ratio_area_s2_modis']=r(tot.s2_ha,tot.modis_ha); out['ratio_co2_s2_modis']=r(tot.s2_tCO2,tot.modis_tCO2)
out['ratio_area_s2_cci']=r(tot.s2_ha,tot.firecci_ha)
out['modis_found_by_s2']=r(tot.s2_and_modis_ha,tot.modis_ha); out['s2_found_by_modis']=r(tot.s2_and_modis_ha,tot.s2_ha)
out['s2_only_ha']=tot.s2_ha-tot.s2_and_modis_ha; out['modis_only_ha']=tot.modis_ha-tot.s2_and_modis_ha
print({k:v for k,v in out.items() if not isinstance(v,dict)})
if h5:
    import h5py
    from matplotlib.path import Path
    aoi=np.array(json.load(open('aoi_coords.json'))[0])
    f=h5py.File(h5[0],'r')
    lat=f['lat'][:]; lon=f['lon'][:]; area=f['ancill/grid_cell_area'][:]
    # fraction of each 0.25 cell inside AOI and north of 33.75, via 10x10 sub-sampling
    sel=(lat[:,0]>33)&(lat[:,0]<39); seln=(lon[0]>-95)&(lon[0]<-86)
    ii=np.where(sel)[0]; jj=np.where(seln)[0]
    p=Path(aoi); frac=np.zeros(lat.shape)
    off=(np.arange(10)+0.5)/10*0.25-0.125
    for i in ii:
        for j in jj:
            la=lat[i,j]+off[:,None]+0*off[None,:]; lo=lon[i,j]+0*off[:,None]+off[None,:]
            pts=np.c_[lo.ravel(),la.ravel()]
            ins=p.contains_points(pts)&(pts[:,1]>=33.75)
            frac[i,j]=ins.mean()
    EF_CO2_AGRI=1585.0  # g CO2 / kg DM, GFED4 agricultural (Akagi et al. 2011)
    g={}
    for m in (9,10,11):
        grp=f['emissions/%02d'%m]
        dm=grp['DM'][:]                       # kg DM m-2 month-1
        agri=grp['partitioning/DM_AGRI'][:]
        dm_agri=dm*agri*area*frac/1e3           # t DM
        dm_all=dm*area*frac/1e3
        g[m]={'DM_agri_t':float(dm_agri.sum()),'CO2_agri_t':float(dm_agri.sum()*EF_CO2_AGRI/1e3),
              'DM_all_t':float(dm_all.sum()),'agri_share':float(dm_agri.sum()/max(dm_all.sum(),1e-9))}
    out['gfed']=g; out['gfed_cells']=int((frac>0).sum())
    print(json.dumps(g,indent=1))
json.dump(out,open('results_subarea_summary.json','w'),indent=1,default=float)

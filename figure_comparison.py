import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
s=json.load(open('results_subarea_summary.json'))
m=['9','10','11']; lab=['September','October','November']
area={'Sentinel-2 + VIIRS (this study)':[s['by_month'][k]['s2_ha'] for k in m],
      'MODIS MCD64A1':[s['by_month'][k]['modis_ha'] for k in m],
      'FireCCI51':[s['by_month'][k]['firecci_ha'] for k in m]}
co2={'Sentinel-2 + VIIRS (this study)':[s['by_month'][k]['s2_tCO2'] for k in m],
     'MODIS MCD64A1':[s['by_month'][k]['modis_tCO2'] for k in m],
     'FireCCI51':[s['by_month'][k]['firecci_tCO2'] for k in m],
     'GFED4.1s agriculture':[s['gfed'][k]['CO2_agri_t'] for k in m]}
col={'Sentinel-2 + VIIRS (this study)':'#d95f0e','MODIS MCD64A1':'#3b6e9c','FireCCI51':'#7fb1d6','GFED4.1s agriculture':'#7a7a7a'}
plt.rcParams.update({'font.size':11})
fig,ax=plt.subplots(1,2,figsize=(12,4.6),dpi=250)
for a,d,yl,t in [(ax[0],area,'Burned cropland (10$^3$ ha)','(a) Burned cropland area'),(ax[1],co2,'CO$_2$ (kt)','(b) CO$_2$ emissions')]:
    n=len(d); w=0.8/n; x=np.arange(3)
    for i,(k,v) in enumerate(d.items()):
        a.bar(x-0.4+w/2+i*w,np.array(v)/1e3,w,label=k,color=col[k],edgecolor='white',lw=0.5)
    a.set_xticks(x); a.set_xticklabels(lab); a.set_ylabel(yl); a.set_title(t,loc='left',fontweight='bold')
    a.spines[['top','right']].set_visible(False); a.grid(axis='y',color='#ddd',lw=0.6); a.set_axisbelow(True)
h,l=ax[1].get_legend_handles_labels()
fig.legend(h,l,loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(0.5,-0.04))
fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig('fig6_comparison.png',bbox_inches='tight')

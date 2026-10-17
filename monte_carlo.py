"""Add Monte Carlo 95% intervals to results.json (same model as the GEE script). usage: mc.py results.json"""
import sys, json, numpy as np
P = {'1':(1.00,1.04,.87,.80),'2':(2.00,2.10,.85,.90),'3':(1.08,1.40,.89,.80),'4':(1.40,1.46,.89,.85),
     '5':(2.10,2.22,.91,.85),'10':(1.00,1.22,.94,.85),'13':(1.00,1.00,.87,.80),'24':(1.30,1.32,.89,.90),
     '28':(1.30,1.46,.89,.85),'29':(1.40,1.50,.90,.85),'45':(0.27,0.47,.85,.80),'46':(0.40,0.42,.22,.85)}
NAME_ID = {'Corn':'1','Cotton':'2','Rice':'3','Sorghum':'4','Soybeans':'5','Peanuts':'10','Popcorn':'13',
           'Winter wheat':'24','Oats':'28','Millet':'29','Sugarcane':'45','Sweet potato':'46'}
d = json.load(open(sys.argv[1]))
ha = {NAME_ID[r.split('|')[0]]: float(r.split('|')[1]) for r in d['crops']}
rng = np.random.default_rng(2020); N = 10000; tot = np.zeros(N); cmc = {}
for k, (lo, hi, dmf, cf) in P.items():
    y = max(float(d['yields'].get(k, 0)), 0)
    a = ha.get(k, 0) * np.maximum(1 + 0.25 * rng.standard_normal(N), 0)
    ef = np.maximum(1 + 0.15 * rng.standard_normal(N), 0) * 1441
    e = a * y * rng.uniform(lo, hi, N) * dmf * cf * ef / 1e6
    cmc[k] = [float(np.percentile(e, 2.5)), float(np.median(e)), float(np.percentile(e, 97.5))]; tot += e
d['crop_mc'] = cmc
d['mc_total'] = [float(np.percentile(tot, 2.5)), float(np.median(tot)), float(np.percentile(tot, 97.5))]
json.dump(d, open(sys.argv[1], 'w'))
print(d['mc_total'])

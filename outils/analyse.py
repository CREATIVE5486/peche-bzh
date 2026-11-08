# Analyse harmonique des marégraphes REFMAR (UTide) + validation + prédiction des PM/BM.
import sys, os, json, glob, datetime as dt
import numpy as np, pandas as pd, utide

OBS, OUT = sys.argv[1], sys.argv[2]
LAT = {"brest":48.383,"lorient":47.7375,"port-tudy":47.644,"crouesty":47.543,"vannes":47.638,
       "larmor-baden":47.58,"le-croisic":47.298,"le-pouliguen":47.27,"saint-nazaire":47.267,"kercabellec":47.43}
SPLIT = pd.Timestamp("2026-09-01")
PRED_FROM, PRED_TO = pd.Timestamp("2026-09-01"), pd.Timestamp("2028-02-01")
HOURLY_FROM, HOURLY_TO = pd.Timestamp("2026-10-01"), pd.Timestamp("2028-01-01")

def load(name):
    df = pd.read_csv(os.path.join(OBS, name+".csv"), names=["t","h"])
    df["t"] = pd.to_datetime(df["t"], format="%Y/%m/%d %H:%M:%S")
    df = df.dropna().drop_duplicates("t").set_index("t").sort_index()
    h = df["h"]
    # suppression des pics aberrants : écart > 0,4 m à la médiane glissante 1 h
    med = h.rolling("60min", center=True).median()
    return h[(h - med).abs() < 0.4]

def extrema(t, h, min_gap_h=2.0):
    """PM/BM par interpolation parabolique sur une série régulière."""
    out = []
    for i in range(1, len(h)-1):
        a, b, c = h[i-1], h[i], h[i+1]
        if not ((b >= a and b > c) or (b <= a and b < c)): continue
        den = a - 2*b + c
        off = 0.5*(a-c)/den if den else 0.0
        step = (t[i+1]-t[i])
        out.append((t[i] + off*step, b - 0.25*(a-c)*off, "PM" if b >= a else "BM"))
    # fusionne les extrema trop proches (tenue du plein / double BM) : garde le plus extrême
    res = []
    for e in out:
        if res and res[-1][2] == e[2] and (e[0]-res[-1][0]) < pd.Timedelta(hours=min_gap_h):
            if (e[2]=="PM" and e[1] > res[-1][1]) or (e[2]=="BM" and e[1] < res[-1][1]): res[-1] = e
            continue
        if res and res[-1][2] != e[2] and (e[0]-res[-1][0]) < pd.Timedelta(hours=min_gap_h):
            continue  # micro-oscillation
        res.append(e)
    return res

def fit(h, lat):
    tt = h.index.to_pydatetime()
    return utide.solve(tt, h.values, lat=lat, method="ols", conf_int="none", trend=False, verbose=False)

def predict(coef, idx):
    return utide.reconstruct(idx.to_pydatetime(), coef, verbose=False).h

def match(obs_e, pred_e, typ, tol_h=3):
    pe = [(t,h) for t,h,k in pred_e if k==typ]
    dts, dhs = [], []
    for t,h,k in obs_e:
        if k != typ: continue
        best = min(pe, key=lambda x: abs((x[0]-t).total_seconds()), default=None)
        if best and abs((best[0]-t).total_seconds()) < tol_h*3600:
            dts.append((best[0]-t).total_seconds()/60); dhs.append(best[1]-h)
    return np.array(dts), np.array(dhs)

report, tides = {}, {}
for name in sorted(LAT):
    p = os.path.join(OBS, name+".csv")
    if not os.path.exists(p) or os.path.getsize(p) < 1000: print(name, "pas de données"); continue
    h = load(name)
    span_days = (h.index[-1]-h.index[0]).days
    train, test = h[h.index < SPLIT], h[h.index >= SPLIT]
    r = {"n": int(len(h)), "debut": str(h.index[0]), "fin": str(h.index[-1]), "jours": span_days,
         "couverture_%": round(100*len(h)/(span_days*144+1),1)}
    if len(train) < 144*200:  # moins de ~200 jours : analyse pas fiable
        r["statut"] = "série trop courte"; report[name] = r; print(name, r); continue
    c = fit(train, LAT[name])
    if len(test) > 144*20:
        tidx = pd.date_range(test.index[0].floor("10min"), test.index[-1], freq="10min")
        ph = pd.Series(predict(c, tidx), index=tidx)
        common = test.index.intersection(tidx)
        resid = test[common] - ph[common]
        r["rms_hauteur_cm"] = round(float(np.sqrt(np.mean((resid - resid.mean())**2)))*100, 1)
        r["biais_cm"] = round(float(resid.mean())*100, 1)
        # extrema observés sur série lissée 30 min (enlève le clapot)
        ts = test.resample("10min").mean().interpolate(limit=3).rolling(3, center=True).mean().dropna()
        oe = extrema(list(ts.index), ts.values)
        pe = extrema(list(tidx), ph.values)
        for typ in ("PM","BM"):
            dtm, dhm = match(oe, pe, typ)
            if len(dtm):
                r[typ] = {"n": int(len(dtm)), "ecart_heure_moy_min": round(float(dtm.mean()),1),
                          "ecart_heure_abs_med_min": round(float(np.median(np.abs(dtm))),1),
                          "ecart_heure_p90_min": round(float(np.percentile(np.abs(dtm),90)),1),
                          "ecart_haut_moy_cm": round(float(dhm.mean())*100,1),
                          "ecart_haut_p90_cm": round(float(np.percentile(np.abs(dhm),90))*100,1)}
    # ajustement final sur toute la série, puis prédiction 2026-2028 au pas de 5 min
    cf = fit(h, LAT[name])
    idx = pd.date_range(PRED_FROM, PRED_TO, freq="5min")
    hp = predict(cf, idx)
    ex = extrema(list(idx), hp)
    r["constituantes"] = len(cf.name)
    r["Z0_m"] = round(float(cf.mean), 3)
    tides[name] = {"z0": float(cf.mean), "ext": [[int(t.timestamp()//60), round(float(v),2), k] for t,v,k in ex if t < HOURLY_TO]}
    sel = (idx >= HOURLY_FROM) & (idx < HOURLY_TO) & (idx.minute == 0)
    tides[name]["hh0"] = int(HOURLY_FROM.timestamp()//60)
    tides[name]["hh"] = [int(round(v*100)) for v in hp[sel]]
    if name == "brest":
        # reconstruction semi-diurne seule (pour le coefficient)
        sd = [n for n in cf.name if n.endswith("2") and n[0] not in "0123456789"]
        hsd = utide.reconstruct(idx.to_pydatetime(), cf, constit=sd, verbose=False).h - cf.mean
        tides[name]["sd_pm"] = [[int(t.timestamp()//60), round(float(v),3)] for t,v,k in extrema(list(idx), hsd) if k=="PM"]
        r["semi_diurnes"] = sd
    report[name] = r
    print(name, json.dumps(r, ensure_ascii=False), flush=True)

json.dump(report, open(os.path.join(OUT,"rapport.json"),"w"), indent=1, ensure_ascii=False)
json.dump(tides, open(os.path.join(OUT,"tides_raw.json"),"w"))

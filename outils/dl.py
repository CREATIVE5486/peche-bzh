# Télécharge les observations REFMAR (1 min) et garde un point toutes les 10 min.
import json, sys, urllib.request, datetime as dt, os, concurrent.futures as cf
OUT = sys.argv[1]
ST = {3:"brest",6349:"lorient",71:"port-tudy",185:"crouesty",347:"vannes",449:"larmor-baden",99:"le-croisic",6421:"le-pouliguen",37:"saint-nazaire",6410:"kercabellec"}
def months():
    d = dt.date(2024,11,1)
    while d < dt.date(2026,11,1):
        n = (d.replace(day=28)+dt.timedelta(days=4)).replace(day=1)
        yield d, n; d = n
def get(sid, a, b):
    url = f"https://services.data.shom.fr/maregraphie/observation/json/{sid}?sources=1&dtStart={a}T00:00:00Z&dtEnd={b}T00:00:00Z"
    for k in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: d = json.load(r)["data"]
            return [(x["timestamp"], x["value"]) for x in d if x["timestamp"][15]=="0" and x["timestamp"][17:]=="00"]
        except Exception as e: err=e
    print("FAIL", sid, a, err, flush=True); return []
def station(sid):
    rows=[]
    for a,b in months(): rows += get(sid, a, b)
    with open(os.path.join(OUT, ST[sid]+".csv"),"w") as f:
        for t,v in rows: f.write(f"{t},{v}\n")
    print(ST[sid], len(rows), flush=True)
with cf.ThreadPoolExecutor(5) as ex: list(ex.map(station, ST))

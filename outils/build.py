# Coefficient (calage sur Brest), décalages port/Brest, Arradon interpolé, contrôle interpolation cosinus, export compact.
import json, sys, numpy as np, datetime as dt
RES, APIM, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
raw = json.load(open(f"{RES}/tides_raw.json")); rep = json.load(open(f"{RES}/rapport.json"))
NAMES = {"brest":"Brest","lorient":"Lorient","port-tudy":"Port-Tudy (Groix)","crouesty":"Le Crouesty",
         "vannes":"Vannes","larmor-baden":"Larmor-Baden","le-croisic":"Le Croisic"}
KEEP = ["lorient","port-tudy","crouesty","vannes","larmor-baden","le-croisic"]

# 1) coefficient : C = A_sd(PM Brest) / U × 100 ; U calé sur les coefs api-maree de Brest (J-30..J+30)
sd = raw["brest"]["sd_pm"]
ref = []
for day in json.load(open(APIM))["data"]:
    for e in day["extrema"]:
        if e["type"]=="PM":
            t = dt.datetime.fromisoformat(day["date"]+"T"+e["time"]).replace(tzinfo=dt.timezone.utc)
            ref.append((int(t.timestamp()//60), e["coef"]))
# PM de Brest complète (pour l'heure), associée au max semi-diurne le plus proche
bpm = [(m,h) for m,h,k in raw["brest"]["ext"] if k=="PM"]
sdt = np.array([m for m,_ in sd]); sda = np.array([a for _,a in sd])
def amp(m): i = np.argmin(abs(sdt-m)); return sda[i] if abs(sdt[i]-m) < 180 else None
pairs = []
for m,c in ref:
    j = min(bpm, key=lambda x: abs(x[0]-m))
    if abs(j[0]-m) <= 30 and amp(j[0]) is not None: pairs.append((amp(j[0]), c))
A = np.array([p[0] for p in pairs]); C = np.array([p[1] for p in pairs])
U = float(np.sum(A*A)/np.sum(A*C))*100  # moindres carrés : C ≈ 100·A/U
err = np.round(100*A/U) - C
print(f"U = {U:.3f} m  sur {len(pairs)} PM ; écart coef moy {np.mean(err):+.2f}, |moy| {np.mean(abs(err)):.2f}, max {np.max(abs(err)):.0f}")
coef = [(m, int(max(20,min(120,round(100*amp(m)/U))))) for m,_ in bpm if amp(m) is not None]

# 2) décalage moyen de la PM locale par rapport à la PM de Brest (min)
def lag(st):
    pm = [m for m,h,k in raw[st]["ext"] if k=="PM"]; bt = np.array([m for m,_ in bpm]); d=[]
    for m in pm:
        i = np.argmin(abs(bt-m)); x = m-bt[i]
        if abs(x) < 240: d.append(x)
    return int(round(np.median(d)))

# 3) contrôle de l'interpolation cosinus entre extrema (vs reconstruction harmonique complète) : fait côté analyse
stations, out = {}, {}
for st in KEEP:
    ext = raw[st]["ext"]
    # alternance stricte PM/BM requise par le format compact
    clean=[ext[0]]
    for e in ext[1:]:
        if e[2]==clean[-1][2]:
            if (e[2]=="PM" and e[1]>clean[-1][1]) or (e[2]=="BM" and e[1]<clean[-1][1]): clean[-1]=e
        else: clean.append(e)
    stations[st] = clean
# Arradon : point milieu Larmor-Baden / Vannes (heure et hauteur)
lb, va = stations["larmor-baden"], stations["vannes"]
arr=[]; vt=np.array([e[0] for e in va])
for e in lb:
    i=np.argmin(abs(vt-e[0])); f=va[i]
    if f[2]==e[2] and abs(f[0]-e[0])<120: arr.append([int(round((e[0]+f[0])/2)), round((e[1]+f[1])/2,2), e[2]])
stations["arradon"]=arr

def hourly(st):
    if st!="arradon": return raw[st]["hh0"], raw[st]["hh"]
    d=(lag("vannes")-lag("larmor-baden"))*60.0
    t0=raw["vannes"]["hh0"]*60; a=np.array(raw["larmor-baden"]["hh"],float); b=np.array(raw["vannes"]["hh"],float)
    tt=t0+3600*np.arange(len(a))
    h=(np.interp(tt-d/2,tt,a)+np.interp(tt+d/2,tt,b))/2
    return raw["vannes"]["hh0"], [int(round(x)) for x in h]
def pack(ext):
    t0=ext[0][0]; prev=t0; dts=[]; hs=[]
    for m,h,k in ext: dts.append(m-prev); prev=m; hs.append(int(round(h*100)))
    dts[0]=0
    return {"t0":t0,"first":ext[0][2],"dt":dts,"h":hs}
data={"stations":{}, "until":"décembre 2027", "coefErr":round(float(np.mean(abs(err))),1), "U":round(U,3)}
for st,ext in stations.items():
    o=pack(ext)
    h0,hh=hourly(st); o["hh0"]=h0; o["hh"]=[hh[0]]+[hh[i]-hh[i-1] for i in range(1,len(hh))]
    if st=="arradon":
        o["lag"]=int(round((lag("larmor-baden")+lag("vannes"))/2))
        r1,r2=rep["larmor-baden"],rep["vannes"]
        o["info"]={"name":"Larmor-Baden + Vannes (moyenne)","months":12,
                   "dtMed":round((r1["PM"]["ecart_heure_abs_med_min"]+r2["PM"]["ecart_heure_abs_med_min"])/2),
                   "dhP90":round((r1["PM"]["ecart_haut_p90_cm"]+r2["PM"]["ecart_haut_p90_cm"])/2)}
    else:
        o["lag"]=lag(st); r=rep[st]
        o["info"]={"name":NAMES[st],"months":round(r["jours"]/30.4),"dtMed":round(r["PM"]["ecart_heure_abs_med_min"]),"dhP90":round(r["PM"]["ecart_haut_p90_cm"])}
    data["stations"][st]=o
    print(st, "décalage PM vs Brest", o["lag"], "min ;", len(ext), "extrema")
c0=coef[0][0]; prev=c0; cd=[]
for m,_ in coef: cd.append(m-prev); prev=m
data["coef"]={"t0":c0,"dt":cd,"c":[c for _,c in coef]}
s=json.dumps(data,separators=(",",":"))
open(OUT,"w").write("window.TIDES_DATA="+s+";")
print("taille", len(s)//1024, "Ko")

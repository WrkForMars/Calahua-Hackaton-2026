import json,base64,io,numpy as np,rasterio,geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.cm as cm
from rasterio.features import rasterize
from indice_infiltracion import CFG,S,N,W,E_,R,nx,ny,T,D
OP="out/parque"; Z=np.load(f"{OP}/subindices_parque.npz")
cls=rasterio.open(f"{OP}/clases_santiago_arteaga.tif").read(1)
raw={"humedad":np.load(f"{D}/humedad.npy"),"ndvi":np.load(f"{D}/ndvi.npy"),"pendiente":Z["slope"],"arena":np.load(f"{D}/arena.npy"),"arcilla":np.load(f"{D}/arcilla.npy"),"fuego":np.load(f"{D}/fuego.npy")}
raw["suelo"]=raw["arena"]-raw["arcilla"]
park=gpd.read_file(f"{OP}/pn_cumbres_conanp.geojson"); mun=gpd.read_file(f"{OP}/municipios_parque.geojson")
foco=mun[mun.municipio.isin(["Santiago","Arteaga"])].copy(); foco["geometry"]=foco.geometry.intersection(park.union_all())
m=(cls>0); rows,cols=np.where(m); r0,r1,c0,c1=rows.min(),rows.max()+1,cols.min(),cols.max()+1
bnd=[[N-r1*R,W+c0*R],[N-r0*R,W+c1*R]]
def png(rgba):
    b=io.BytesIO(); plt.imsave(b,rgba[r0:r1,c0:c1],format="png"); return "data:image/png;base64,"+base64.b64encode(b.getvalue()).decode()
def cont(a,cmap,lo,hi):
    rgba=(matplotlib.colormaps[cmap](np.clip((a-lo)/(hi-lo),0,1))*255).astype(np.uint8); rgba[...,3]=np.where(m&np.isfinite(a),210,0); return rgba
def ramp(cmap): return "linear-gradient(to right,"+",".join(matplotlib.colors.to_hex(matplotlib.colormaps[cmap](t)) for t in np.linspace(0,1,6))+")"
crgba=np.zeros((ny,nx,4),np.uint8)
for c,v in {1:(255,237,160),2:(253,141,60),3:(189,0,38)}.items(): crgba[cls==c]=(*v,220)
layers=[{"id":"clases","nombre":"Prioridad (alta/media/baja)","img":png(crgba),"show":True}]
from rasterio.warp import reproject,Resampling
Gd=rasterio.open("../resultados/clases_gee_v2.tif"); gee=np.zeros_like(cls)
reproject(Gd.read(1),gee,src_transform=Gd.transform,src_crs=Gd.crs,dst_transform=T,dst_crs="EPSG:4326",resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
grgba=np.zeros((ny,nx,4),np.uint8)
for c,v in {1:(255,237,160),2:(253,141,60),3:(189,0,38)}.items(): grgba[gee==c]=(*v,220)
_m=m.copy(); m=m|(gee>0); r_=np.where(m); r0,r1,c0,c1=r_[0].min(),r_[0].max()+1,r_[1].min(),r_[1].max()+1; bnd=[[N-r1*R,W+c0*R],[N-r0*R,W+c1*R]]
layers[0]["img"]=png(crgba)
layers.append({"id":"gee","nombre":"Validación: clases GEE v2","img":png(grgba),"show":False})
m=_m
spec=[("humedad","Humedad del suelo (TerraClimate)","mm","YlGnBu"),("ndvi","Vegetación – NDVI (MODIS 2023–24)","adim.","YlGn"),("pendiente","Pendiente (Copernicus DEM)","°","magma_r"),("suelo","Suelo: arena − arcilla (SoilGrids 0–30 cm)","puntos %","BrBG_r"),("fuego","Incendios: meses quemados (MCD64A1 2012–2025)","meses","OrRd")]
for k,nm,u,cmap in spec:
    v=raw[k][m]; lo,hi=(0,max(1,float(np.nanmax(v)))) if k=="fuego" else tuple(float(x) for x in np.nanpercentile(v,[2,98]))
    layers.append({"id":k,"nombre":nm,"unidad":u,"img":png(cont(raw[k],cmap,lo,hi)),"lo":round(lo,2),"hi":round(hi,2),"ramp":ramp(cmap),"show":False})
# grid de consulta (0.002°)
st=2; sub=lambda a: a[r0:r1:st,c0:c1:st]
keys={"hum":raw["humedad"],"ndvi":raw["ndvi"],"pend":raw["pendiente"],"arena":raw["arena"],"arcilla":raw["arcilla"],"fuego":raw["fuego"],
      "s_hum":Z["s_humedad"],"s_veg":Z["s_vegetacion"],"s_pen":Z["s_pendiente"],"s_sue":Z["s_suelo"],"s_inc":Z["s_incendios"],"idx":Z["idx"]}
mm=sub(m); grid={"lat0":N-r0*R,"lon0":W+c0*R,"res":R*st,"nr":int(mm.shape[0]),"nc":int(mm.shape[1]),"cls":sub(cls).astype(int).tolist()}
for k,a in keys.items(): grid[k]=[[None if (not mm[i,j] or not np.isfinite(x)) else round(float(x),3 if k in("ndvi","idx") or k.startswith("s_") else 1) for j,x in enumerate(row)] for i,row in enumerate(sub(a))]
# resúmenes
pix=(R*111.32)**2*np.cos(np.radians((S+N)/2)); munid_r=rasterize([(g,i+1) for i,g in enumerate(foco.geometry)],out_shape=(ny,nx),transform=T)
cards={}
for i,n in enumerate(foco.municipio,1):
    reg=m&(munid_r==i); c={"km2":round(float(reg.sum()*pix),1),"clases":{}}
    for cv,cn in [(3,"alta"),(2,"media"),(1,"baja")]:
        rr=reg&(cls==cv); d={"km2":round(float(rr.sum()*pix),1),"pct":round(100*float(rr.sum()/reg.sum()),1)}
        if rr.sum():
            for k in ["humedad","ndvi","pendiente","arena","arcilla","fuego"]: d[k]=round(float(np.nanmean(raw[k][rr])),3 if k=="ndvi" else 2)
            d["indice"]=round(float(np.nanmean(Z["idx"][rr])),3)
        c["clases"][cn]=d
    cards[n]=c
json.dump(cards,open(f"{OP}/resumen_pagina.json","w"),indent=1,ensure_ascii=False)
W_=CFG["weights"]
cmp_=json.load(open("out/parque/comparacion_gee_v2.json")); st_=json.load(open("out/parque/estadisticas_parque.json"))
VAL={"cmp":cmp_,"py_terciles":[round(x,4) for x in st_["santiago_arteaga"]["umbrales"]],"gee_terciles":[0.36519,0.50585],
     "gee_km2_consola":{"baja":176.31,"media":177.65,"alta":178.24}}
html=open("plantilla_pagina.html",encoding="utf-8").read()
for k,v in {"__LAYERS__":json.dumps(layers),"__BOUNDS__":json.dumps(bnd),"__GRID__":json.dumps(grid,separators=(",",":")),"__CARDS__":json.dumps(cards,ensure_ascii=False),
            "__PARK__":park.to_json(),"__MUN__":foco.to_json(),"__W__":json.dumps(W_),"__VAL__":json.dumps(VAL,ensure_ascii=False)}.items(): html=html.replace(k,v)
open(f"{OP}/pagina_santiago_arteaga.html","w",encoding="utf-8").write(html); print(json.dumps(cards,indent=1,ensure_ascii=False))

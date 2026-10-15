"""Índice de prioridad recortado al Parque Nacional Cumbres de Monterrey (polígono oficial CONANP),
con énfasis en Santiago y Arteaga. Reutiliza capas/caché de indice_infiltracion.py."""
import json, numpy as np, rasterio, geopandas as gpd
from rasterio.features import rasterize
import indice_infiltracion as B
from indice_infiltracion import CFG, W, S, E_, N, R, nx, ny, T, CRS, D, O
import os; OP=f"{O}/parque"; os.makedirs(OP,exist_ok=True)
PARK_F=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","datos","pn_cumbres.geojson")
park=gpd.read_file(PARK_F).to_crs(4326)
if not os.path.exists(f"{D}/mex_adm2.geojson"):  # geoBoundaries MEX ADM2 (CC BY 4.0)
    import requests; open(f"{D}/mex_adm2.geojson","wb").write(requests.get("https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/MEX/ADM2/geoBoundaries-MEX-ADM2.geojson",timeout=600).content)
mun=gpd.read_file(f"{D}/mex_adm2.geojson").to_crs(4326)
mun=mun[mun.intersects(park.union_all())].rename(columns={"shapeName":"municipio"})[["municipio","geometry"]]
FOCO=CFG.get("municipios_foco",["Santiago","Arteaga"])
inpark=rasterize(park.geometry,out_shape=(ny,nx),transform=T,fill=0,default_value=1,all_touched=False).astype(bool)
munid=rasterize([(g,i+1) for i,g in enumerate(mun.geometry)],out_shape=(ny,nx),transform=T,fill=0).astype(np.int16)
foco_ids=[i+1 for i,n in enumerate(mun.municipio) if n in FOCO]
infoco=inpark&np.isin(munid,foco_ids)

L={"dem":B.cache("dem",B.dem)}
for k,fn in [("ndvi",B.ndvi),("fuego",B.fuego),("arena",B.soil("sand")),("arcilla",B.soil("clay")),("humedad",B.humedad)]:
    try: L[k]=B.cache(k,fn)
    except Exception as e: print("FALLÓ",k,e)
z=L["dem"]; dy=R*111320; dx=R*111320*np.cos(np.radians((S+N)/2))
gy,gx=np.gradient(z,dy,dx); slope=np.degrees(np.arctan(np.hypot(gx,gy)))
def nrm(a,reg):  # normalización con percentiles 2–98 calculados DENTRO del parque
    lo,hi=np.nanpercentile(a[reg],[2,98]); return np.clip((a-lo)/(hi-lo),0,1)
sub={"pendiente":np.where(slope<=CFG["slope_opt_max"],1.0,np.clip(1-(slope-CFG["slope_opt_max"])/(CFG["slope_max"]-CFG["slope_opt_max"]),0,1))}
if "ndvi" in L: sub["vegetacion"]=1-nrm(L["ndvi"],inpark)
if "humedad" in L: sub["humedad"]=1-nrm(L["humedad"],inpark)
if "fuego" in L: sub["incendios"]=np.clip(L["fuego"]/max(1,L["fuego"][inpark].max()),0,1)
if "arena" in L and "arcilla" in L: sub["suelo"]=nrm(L["arena"]-L["arcilla"],inpark)
w={k:v for k,v in CFG["weights"].items() if k in sub}; tot=sum(w.values())
idx=sum(np.nan_to_num(sub[k],nan=np.nan)*w[k] for k in w)/tot
valid=inpark&np.isfinite(idx); idx=np.where(valid,idx,np.nan)
def classify(reg):
    t1,t2=np.nanpercentile(idx[reg],[33.333,66.667])
    return np.where(~reg,0,np.where(idx>=t2,3,np.where(idx>=t1,2,1))).astype(np.uint8),(float(t1),float(t2))
np.savez(f"{OP}/subindices_parque.npz",idx=idx,slope=slope,infoco=infoco,munid=munid,**{"s_"+k:v for k,v in sub.items()})
cls_p,tp=classify(valid); cls_f,tf=classify(valid&infoco)
pix=(R*111.32)**2*np.cos(np.radians((S+N)/2))
def areas(cls,reg):
    n=reg.sum(); return {nm:{"km2":round(float((cls==c).sum()*pix),1),"pct":round(100*float((cls==c).sum())/n,1)} for c,nm in [(3,"alta"),(2,"media"),(1,"baja")]}
st={"poligono":"CONANP shp_anp_20260924 – PN Cumbres de Monterrey (NUM_ANP 144, decreto DOF 2000-11-17, 177,395.95 ha)",
    "capas_usadas":list(w),"pesos_normalizados":{k:round(v/tot,3) for k,v in w.items()},
    "parque":{"km2_rasterizado":round(float(valid.sum()*pix),1),"umbrales":tp,"clases":areas(cls_p,valid)},
    "santiago_arteaga":{"km2_rasterizado":round(float((valid&infoco).sum()*pix),1),"umbrales":tf,"clases":areas(cls_f,valid&infoco)},
    "por_municipio_(terciles_del_parque)":{}}
for i,n in enumerate(mun.municipio,1):
    r=valid&(munid==i)
    if r.sum()*pix>1: st["por_municipio_(terciles_del_parque)"][n]={"km2":round(float(r.sum()*pix),1),"indice_medio":round(float(np.nanmean(idx[r])),3),"km2_alta":round(float(((cls_p==3)&r).sum()*pix),1),"pct_alta":round(100*float(((cls_p==3)&r).sum()/r.sum()),1)}
# subzonas: celdas de ~5 km con mayor % de clase alta (terciles Santiago+Arteaga) + localidad más cercana (lista en config)
locs=CFG["localidades_parque"]; cells=[]
for i in range(0,ny,50):
    for j in range(0,nx,50):
        reg=(valid&infoco)[i:i+50,j:j+50]
        if reg.sum()>800:
            lo=W+(j+25)*R; la=N-(i+25)*R
            near=min(locs,key=lambda k:(locs[k][0]-lo)**2+(locs[k][1]-la)**2)
            dkm=(((locs[near][0]-lo)*101)**2+((locs[near][1]-la)*111)**2)**.5
            cells.append({"lon":round(lo,3),"lat":round(la,3),"pct_alta":round(100*float((cls_f[i:i+50,j:j+50]==3).sum()/reg.sum()),1),"indice_medio":round(float(np.nanmean(idx[i:i+50,j:j+50][reg])),3),"elev_media_m":int(np.nanmean(z[i:i+50,j:j+50][reg])),"localidad_cercana":f"{near} (~{dkm:.0f} km)"})
st["top_subzonas_santiago_arteaga_(celdas_5km)"]=sorted(cells,key=lambda c:-c["pct_alta"])[:10]
json.dump(st,open(f"{OP}/estadisticas_parque.json","w"),indent=2,ensure_ascii=False); print(json.dumps(st,indent=1,ensure_ascii=False))
prof=dict(driver="GTiff",width=nx,height=ny,crs=CRS,transform=T,compress="deflate")
cmap={1:(255,237,160),2:(253,141,60),3:(189,0,38)}
with rasterio.open(f"{OP}/indice_parque.tif","w",count=1,dtype="float32",nodata=np.nan,**prof) as d: d.write(idx.astype(np.float32),1)
for nm,c in [("clases_parque.tif",cls_p),("clases_santiago_arteaga.tif",cls_f)]:
    with rasterio.open(f"{OP}/{nm}","w",count=1,dtype="uint8",nodata=0,**prof) as d: d.write(c,1); d.write_colormap(1,cmap)
park.to_file(f"{OP}/pn_cumbres_conanp.geojson"); mun.to_file(f"{OP}/municipios_parque.geojson")
# PNG
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LightSource; from matplotlib.patches import Patch; from matplotlib.lines import Line2D
b=park.total_bounds; pad=.05; ext=[W,E_,S,N]
hs=LightSource(315,45).hillshade(np.nan_to_num(z),vert_exag=.05,dx=dx,dy=dy)
fig,axs=plt.subplots(1,2,figsize=(20,9.5))
for ax,cls,ttl,ar in [(axs[0],cls_p,"Terciles dentro de todo el parque",st["parque"]),(axs[1],cls_f,"Terciles solo en Santiago + Arteaga (dentro del parque)",st["santiago_arteaga"])]:
    ax.imshow(hs,cmap="gray",extent=ext)
    ax.imshow(np.ma.masked_equal(cls,0),cmap=ListedColormap(["#ffeda0","#fd8d3c","#bd0026"]),vmin=1,vmax=3,alpha=.8,extent=ext,interpolation="nearest")
    mun.boundary.plot(ax=ax,color="#333",lw=.6,ls="--")
    mun[mun.municipio.isin(FOCO)].boundary.plot(ax=ax,color="#0047ab",lw=2)
    park.boundary.plot(ax=ax,color="#006400",lw=2)
    for _,r in mun.iterrows():
        c=r.geometry.intersection(park.union_all())
        if not c.is_empty and c.area>1e-4: p=c.representative_point(); ax.text(p.x,p.y,r.municipio,fontsize=9,ha="center",color="#0047ab" if r.municipio in FOCO else "#222",weight="bold" if r.municipio in FOCO else "normal",bbox=dict(fc="white",alpha=.6,lw=0))
    ax.set_xlim(b[0]-pad,b[2]+pad); ax.set_ylim(b[1]-pad,b[3]+pad)
    ax.set_title(ttl,fontsize=12); ax.set_xlabel("Longitud"); ax.set_ylabel("Latitud")
    ax.legend(handles=[Patch(color="#bd0026",label=f"Alta ({ar['clases']['alta']['km2']} km²)"),Patch(color="#fd8d3c",label=f"Media ({ar['clases']['media']['km2']} km²)"),Patch(color="#ffeda0",label=f"Baja ({ar['clases']['baja']['km2']} km²)"),
        Line2D([],[],color="#006400",lw=2,label="PN Cumbres de Monterrey (CONANP)"),Line2D([],[],color="#0047ab",lw=2,label="Santiago y Arteaga"),Line2D([],[],color="#333",lw=.6,ls="--",label="Otros municipios")],loc="lower right",fontsize=9,title="Prioridad de restauración")
fig.suptitle("Reto 1 · Infiltración — Prioridad de restauración de suelos en el Parque Nacional Cumbres de Monterrey",fontsize=15)
fig.text(.5,.01,"Pesos: "+", ".join(f"{k} {v:.0%}" for k,v in st["pesos_normalizados"].items())+" · Fuentes: CONANP (ANP 2026), geoBoundaries/INEGI (municipios 2012), Copernicus DEM, MODIS MOD13Q1/MCD64A1, SoilGrids 2.0, TerraClimate",ha="center",fontsize=9)
plt.tight_layout(rect=[0,.03,1,.95]); plt.savefig(f"{OP}/mapa_parque.png",dpi=140)
# HTML
import folium, base64, io
def ov(cls):
    rgba=np.zeros((ny,nx,4),np.uint8)
    for c,v in cmap.items(): rgba[cls==c]=(*v,200)
    buf=io.BytesIO(); plt.imsave(buf,rgba,format="png"); return "data:image/png;base64,"+base64.b64encode(buf.getvalue()).decode()
m=folium.Map(location=[(b[1]+b[3])/2,(b[0]+b[2])/2],zoom_start=10)
folium.TileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",attr="Esri",name="Satélite").add_to(m)
folium.raster_layers.ImageOverlay(ov(cls_p),bounds=[[S,W],[N,E_]],name="Prioridad (terciles parque)").add_to(m)
folium.raster_layers.ImageOverlay(ov(cls_f),bounds=[[S,W],[N,E_]],name="Prioridad (terciles Santiago+Arteaga)",show=False).add_to(m)
folium.GeoJson(park,name="PN Cumbres de Monterrey",style_function=lambda f:{"color":"#006400","weight":3,"fill":False}).add_to(m)
folium.GeoJson(mun,name="Municipios",style_function=lambda f:{"color":"#0047ab" if f["properties"]["municipio"] in FOCO else "#444","weight":3 if f["properties"]["municipio"] in FOCO else 1,"dashArray":None if f["properties"]["municipio"] in FOCO else "5,5","fill":False},tooltip=folium.GeoJsonTooltip(["municipio"])).add_to(m)
for c in st["top_subzonas_santiago_arteaga_(celdas_5km)"][:5]:
    folium.CircleMarker([c["lat"],c["lon"]],radius=10,color="black",fill=False,tooltip=f"{c['pct_alta']}% alta · cerca de {c['localidad_cercana']} · {c['elev_media_m']} m").add_to(m)
folium.LayerControl(collapsed=False).add_to(m)
m.get_root().html.add_child(folium.Element('<div style="position:fixed;bottom:30px;left:30px;z-index:9999;background:white;padding:8px;border:1px solid #888;font-size:13px"><b>Prioridad de restauración</b><br><span style="background:#bd0026">&nbsp;&nbsp;&nbsp;</span> Alta<br><span style="background:#fd8d3c">&nbsp;&nbsp;&nbsp;</span> Media<br><span style="background:#ffeda0">&nbsp;&nbsp;&nbsp;</span> Baja<br><span style="color:#006400">━</span> PN Cumbres (CONANP)<br><span style="color:#0047ab">━</span> Santiago / Arteaga<br>○ Top 5 subzonas</div>'))
m.save(f"{OP}/mapa_parque.html"); print("Listo.")

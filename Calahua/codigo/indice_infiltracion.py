"""Índice de prioridad de restauración de suelos para infiltración – partes altas ZMM (Lumi)."""
import json, os, numpy as np, rasterio, warnings
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from concurrent.futures import ThreadPoolExecutor
import pystac_client, planetary_computer as pc
warnings.filterwarnings("ignore")
BASE=os.path.dirname(os.path.abspath(__file__)); D=f"{BASE}/data"; O=f"{BASE}/out"
os.makedirs(D,exist_ok=True); os.makedirs(O,exist_ok=True)
CFG=json.load(open(f"{BASE}/config.json"))
W,S,E_,N=CFG["bbox"]; R=CFG["res_deg"]
nx,ny=round((E_-W)/R),round((N-S)/R); T=from_origin(W,N,R,R); CRS="EPSG:4326"
def warp(src_path_or_ds, band=1, rs=Resampling.bilinear, nodata=None):
    out=np.full((ny,nx),np.nan,np.float32)
    ds=rasterio.open(src_path_or_ds) if isinstance(src_path_or_ds,str) else src_path_or_ds
    reproject(rasterio.band(ds,band),out,dst_transform=T,dst_crs=CRS,resampling=rs,
              src_nodata=nodata if nodata is not None else ds.nodata,dst_nodata=np.nan)
    return out
def cache(name,fn):
    p=f"{D}/{name}.npy"
    if os.path.exists(p): return np.load(p)
    a=fn(); np.save(p,a); return a
STAC=pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=pc.sign_inplace)

# 1) DEM Copernicus GLO-30 (AWS open data)
def dem():
    acc=np.full((ny,nx),np.nan,np.float32)
    for lat in range(int(np.floor(S)),int(np.ceil(N))):
        for lon in range(int(np.floor(W)),int(np.ceil(E_))):
            t=f"Copernicus_DSM_COG_10_N{lat:02d}_00_W{abs(lon):03d}_00_DEM"
            a=warp(f"https://copernicus-dem-30m.s3.amazonaws.com/{t}/{t}.tif",rs=Resampling.average)
            acc=np.where(np.isnan(acc),a,acc)
    return acc
# 2) NDVI MODIS MOD13Q1 v061 (mediana de compuestos 16 días)
def ndvi():
    its=list(STAC.search(collections=["modis-13Q1-061"],bbox=CFG["bbox"],datetime=CFG["ndvi_period"]).items())
    its=[i for i in its if i.id.startswith("MOD13Q1")]
    def rd(i):
        a=warp(i.assets["250m_16_days_NDVI"].href,rs=Resampling.average,nodata=-3000)*1e-4
        return i.properties["start_datetime"][:10],a
    with ThreadPoolExecutor(8) as ex: res=list(ex.map(rd,its))
    by={}
    for d,a in res: by.setdefault(d,np.full((ny,nx),np.nan,np.float32)); by[d]=np.where(np.isnan(by[d]),a,by[d])
    print(f"  NDVI: {len(by)} compuestos, {len(its)} items")
    return np.nanmedian(np.stack(list(by.values())),0)
# 3) Área quemada MODIS MCD64A1 v061 (nº de meses con quema)
def fuego():
    its=list(STAC.search(collections=["modis-64A1-061"],bbox=CFG["bbox"],datetime=CFG["fire_period"]).items())
    def rd(i):
        a=warp(i.assets["Burn_Date"].href,rs=Resampling.nearest,nodata=-9999)
        return (a>0).astype(np.float32)
    with ThreadPoolExecutor(8) as ex: res=list(ex.map(rd,its))
    print(f"  Fuego: {len(its)} items mensuales")
    return np.sum(res,0)
# 4) Suelo SoilGrids 2.0 (ISRIC WCS) arena y arcilla 0-30 cm
def soil(prop):
    def f():
        layers=[]
        for dep in ["0-5cm","5-15cm","15-30cm"]:
            p=f"{D}/{prop}_{dep}.tif"
            if not os.path.exists(p):
                u=(f"https://maps.isric.org/mapserv?map=/map/{prop}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage"
                   f"&COVERAGEID={prop}_{dep}_mean&FORMAT=GEOTIFF_INT16&SUBSET=X({W},{E_})&SUBSET=Y({S},{N})"
                   "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
                import requests; r=requests.get(u,timeout=300); r.raise_for_status(); open(p,"wb").write(r.content)
            a=warp(p,nodata=-32768); a[a<=0]=np.nan; layers.append(a/10)  # g/kg -> %
        return np.nanmean(layers,0)
    return f
# 5) Humedad del suelo: TerraClimate (mm, columna de suelo) promedio
def humedad():
    import xarray as xr
    col=STAC.get_collection("terraclimate"); a=col.assets["zarr-abfs"]
    ds=xr.open_dataset(a.href,**a.extra_fields["xarray:open_kwargs"])
    y0,y1=CFG["sm_years"]
    v=ds["soil"].sel(lon=slice(W-.1,E_+.1),lat=slice(N+.1,S-.1),time=slice(f"{y0}",f"{y1}")).mean("time").load()
    lat=v.lat.values; lon=v.lon.values; arr=v.values.astype(np.float32)
    src_T=from_origin(lon[0]-(lon[1]-lon[0])/2,lat[0]+(lat[0]-lat[1])/2,lon[1]-lon[0],lat[0]-lat[1])
    out=np.full((ny,nx),np.nan,np.float32)
    reproject(arr,out,src_transform=src_T,src_crs=CRS,dst_transform=T,dst_crs=CRS,resampling=Resampling.bilinear,src_nodata=np.nan,dst_nodata=np.nan)
    return out

def norm(a,lo=None,hi=None):
    lo=np.nanpercentile(a,2) if lo is None else lo; hi=np.nanpercentile(a,98) if hi is None else hi
    return np.clip((a-lo)/(hi-lo),0,1)

if __name__=="__main__":
    print("Descargando/cargando capas…")
    layers={}; z=cache("dem",dem); layers["dem"]=z
    for k,fn in [("ndvi",ndvi),("fuego",fuego),("arena",soil("sand")),("arcilla",soil("clay")),("humedad",humedad)]:
        try: layers[k]=cache(k,fn); print(f"  {k}: OK")
        except Exception as e: print(f"  {k}: FALLÓ -> {e}")
    # pendiente (grados)
    dy=R*111320; dx=R*111320*np.cos(np.radians((S+N)/2))
    gy,gx=np.gradient(z,dy,dx); slope=np.degrees(np.arctan(np.hypot(gx,gy)))
    # sub-índices 0-1 (1 = mayor prioridad)
    sub={}
    sub["pendiente"]=np.where(slope<=CFG["slope_opt_max"],1.0,np.clip(1-(slope-CFG["slope_opt_max"])/(CFG["slope_max"]-CFG["slope_opt_max"]),0,1))
    if "ndvi" in layers: sub["vegetacion"]=1-norm(layers["ndvi"])          # vegetación degradada
    if "humedad" in layers: sub["humedad"]=1-norm(layers["humedad"])     # déficit de humedad
    if "fuego" in layers: sub["incendios"]=norm(layers["fuego"],0,max(1,np.nanpercentile(layers["fuego"],99.5)))
    if "arena" in layers and "arcilla" in layers: sub["suelo"]=norm(layers["arena"]-layers["arcilla"])  # mayor capacidad de infiltración
    w={k:v for k,v in CFG["weights"].items() if k in sub}; tot=sum(w.values())
    idx=sum(sub[k]*w[k] for k in w)/tot
    mask=(z>=CFG["min_elev_m"])&np.isfinite(idx)
    idx=np.where(mask,idx,np.nan)
    t1,t2=np.nanpercentile(idx,[33.333,66.667])
    cls=np.where(~mask,0,np.where(idx>=t2,3,np.where(idx>=t1,2,1))).astype(np.uint8)
    prof=dict(driver="GTiff",width=nx,height=ny,crs=CRS,transform=T,compress="deflate")
    with rasterio.open(f"{O}/indice_prioridad.tif","w",count=1,dtype="float32",nodata=np.nan,**prof) as d: d.write(idx.astype(np.float32),1); d.set_band_description(1,"indice 0-1")
    with rasterio.open(f"{O}/clases_prioridad.tif","w",count=1,dtype="uint8",nodata=0,**prof) as d: d.write(cls,1); d.write_colormap(1,{1:(255,237,160),2:(253,141,60),3:(189,0,38)})
    with rasterio.open(f"{O}/subindices.tif","w",count=len(sub),dtype="float32",nodata=np.nan,**prof) as d:
        for i,(k,v) in enumerate(sub.items(),1): d.write(np.where(mask,v,np.nan).astype(np.float32),i); d.set_band_description(i,k)
    # estadísticas
    pix_km2=(R*111.32)**2*np.cos(np.radians((S+N)/2)); n=mask.sum()
    stats={"capas_usadas":list(w),"pesos_normalizados":{k:round(v/tot,3) for k,v in w.items()},
           "umbrales_terciles":[float(t1),float(t2)],"area_analizada_km2":round(float(n*pix_km2),1),
           "clases":{nm:{"km2":round(float((cls==c).sum()*pix_km2),1),"pct":round(100*float((cls==c).sum())/n,1)} for c,nm in [(3,"alta"),(2,"media"),(1,"baja")]}}
    zones={}
    for nm,(lo,la) in CFG["zonas_referencia"].items():
        r=int((N-la)/R); c=int((lo-W)/R); win=(slice(max(r-50,0),r+50),slice(max(c-50,0),c+50))
        m=mask[win]
        if m.sum()>0: zones[nm]={"indice_medio":round(float(np.nanmean(idx[win])),3),"pct_alta":round(100*float((cls[win]==3).sum()/m.sum()),1),"km2_analizados":round(float(m.sum()*pix_km2),1)}
    stats["zonas_referencia_(~10x10km)"]=dict(sorted(zones.items(),key=lambda kv:-kv[1]["pct_alta"]))
    for k,a in layers.items(): stats.setdefault("rangos_capas",{})[k]=[round(float(np.nanmin(a)),3),round(float(np.nanmax(a)),3)]
    json.dump(stats,open(f"{O}/estadisticas.json","w"),indent=2,ensure_ascii=False); print(json.dumps(stats,indent=1,ensure_ascii=False))
    # mapa PNG
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, LightSource; from matplotlib.patches import Patch
    hs=LightSource(315,45).hillshade(np.nan_to_num(z),vert_exag=0.05,dx=dx,dy=dy)
    fig,ax=plt.subplots(figsize=(12,10.5)); ext=[W,E_,S,N]
    ax.imshow(hs,cmap="gray",extent=ext)
    ax.imshow(np.ma.masked_equal(cls,0),cmap=ListedColormap(["#ffeda0","#fd8d3c","#bd0026"]),vmin=1,vmax=3,alpha=.7,extent=ext,interpolation="nearest")
    for nm,(lo,la) in CFG["lugares"].items(): ax.plot(lo,la,"k^",ms=6); ax.annotate(nm,(lo,la),xytext=(4,4),textcoords="offset points",fontsize=9,weight="bold")
    ax.legend(handles=[Patch(color="#bd0026",label=f"Alta ({stats['clases']['alta']['pct']}%)"),Patch(color="#fd8d3c",label=f"Media ({stats['clases']['media']['pct']}%)"),Patch(color="#ffeda0",label=f"Baja ({stats['clases']['baja']['pct']}%)"),Patch(color="#999",label=f"Fuera de análisis (< {CFG['min_elev_m']} m)")],loc="lower left",title="Prioridad de restauración",fontsize=10)
    ax.set_title("Reto 1 · Infiltración — ¿Dónde restaurar suelos en partes altas?\nÍndice de prioridad (terciles) · Sierra Madre Oriental, cuencas altas ZMM",fontsize=13)
    ax.set_xlabel("Longitud"); ax.set_ylabel("Latitud")
    ax.text(.99,.01,"Pesos: "+", ".join(f"{k} {v:.0%}" for k,v in stats["pesos_normalizados"].items())+"\nFuentes: Copernicus DEM GLO-30, MODIS MOD13Q1/MCD64A1, SoilGrids 2.0, TerraClimate",transform=ax.transAxes,ha="right",va="bottom",fontsize=8,bbox=dict(fc="white",alpha=.8))
    ax.set_xlim(W,E_); ax.set_ylim(S,N); plt.tight_layout(); plt.savefig(f"{O}/mapa_prioridad.png",dpi=150)
    # HTML folium
    import folium, base64, io
    rgba=np.zeros((ny,nx,4),np.uint8); cols={1:(255,237,160),2:(253,141,60),3:(189,0,38)}
    for c,v in cols.items(): rgba[cls==c]=(*v,190)
    buf=io.BytesIO(); plt.imsave(buf,rgba,format="png")
    m=folium.Map(location=[(S+N)/2,(W+E_)/2],zoom_start=9,tiles="OpenStreetMap")
    folium.TileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",attr="Esri",name="Satélite").add_to(m)
    folium.raster_layers.ImageOverlay("data:image/png;base64,"+base64.b64encode(buf.getvalue()).decode(),bounds=[[S,W],[N,E_]],name="Prioridad (alta/media/baja)",opacity=.75).add_to(m)
    for nm,(lo,la) in CFG["lugares"].items(): folium.Marker([la,lo],tooltip=nm).add_to(m)
    folium.LayerControl().add_to(m)
    m.get_root().html.add_child(folium.Element('<div style="position:fixed;bottom:30px;left:30px;z-index:9999;background:white;padding:8px;border:1px solid #888;font-size:13px"><b>Prioridad de restauración</b><br><span style="background:#bd0026">&nbsp;&nbsp;&nbsp;</span> Alta<br><span style="background:#fd8d3c">&nbsp;&nbsp;&nbsp;</span> Media<br><span style="background:#ffeda0">&nbsp;&nbsp;&nbsp;</span> Baja</div>'))
    m.save(f"{O}/mapa_prioridad.html"); print("Listo.")

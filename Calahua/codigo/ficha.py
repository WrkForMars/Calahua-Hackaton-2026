import json,numpy as np,rasterio,geopandas as gpd
from rasterio.features import rasterize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap; from matplotlib.patches import Patch
from indice_infiltracion import S,N,W,E_,R,nx,ny,T,D
OP="out/parque"
cls=rasterio.open(f"{OP}/clases_santiago_arteaga.tif").read(1); idx=rasterio.open(f"{OP}/indice_parque.tif").read(1)
z=np.load(f"{D}/dem.npy"); ndvi=np.load(f"{D}/ndvi.npy")
dy=R*111320; dx=R*111320*np.cos(np.radians((S+N)/2)); gy,gx=np.gradient(z,dy,dx); sl=np.degrees(np.arctan(np.hypot(gx,gy)))
park=gpd.read_file(f"{OP}/pn_cumbres_conanp.geojson"); mun=gpd.read_file(f"{OP}/municipios_parque.geojson")
pix=(R*111.32)**2*np.cos(np.radians((S+N)/2)); st={}
fig=plt.figure(figsize=(8.27,11.69)); ext=[W,E_,S,N]
fig.text(.5,.965,"Ficha · Reto 1 Infiltración",ha="center",fontsize=16,weight="bold")
fig.text(.5,.945,"Prioridad de restauración de suelos — PN Cumbres de Monterrey: porciones de Santiago y Arteaga",ha="center",fontsize=10.5)
for k,(name,y0) in enumerate([("Santiago",.555),("Arteaga",.17)]):
    g=mun[mun.municipio==name].geometry.intersection(park.union_all()).union_all()
    m=rasterize([g],out_shape=(ny,nx),transform=T).astype(bool)&(cls>0)
    n=m.sum(); d={"km2":n*pix}
    for c,nm in [(3,"alta"),(2,"media"),(1,"baja")]: d[nm]=((cls==c)&m).sum()*pix; d[nm+"_pct"]=100*((cls==c)&m).sum()/n
    d.update(indice=np.nanmean(idx[m]),elev=np.nanmean(z[m]),pend=np.nanmean(sl[m]),ndvi=np.nanmean(ndvi[m])); st[name]={a:round(float(b),3) for a,b in d.items()}
    ax=fig.add_axes([.07,y0,.46,.34]); ax.set_anchor("N"); b=g.bounds; p=.02
    show=np.where(m,cls,0); ax.imshow(np.ma.masked_equal(show,0),cmap=ListedColormap(["#ffeda0","#fd8d3c","#bd0026"]),vmin=1,vmax=3,extent=ext,interpolation="nearest")
    gpd.GeoSeries([g]).boundary.plot(ax=ax,color="#0047ab",lw=1.2); park.boundary.plot(ax=ax,color="#006400",lw=.8,ls="--")
    ax.set_xlim(b[0]-p,b[2]+p); ax.set_ylim(b[1]-p,b[3]+p); ax.tick_params(labelsize=7); ax.set_xlabel("Longitud",fontsize=8); ax.set_ylabel("Latitud",fontsize=8); ax.set_title(f"{name} (dentro del parque)",fontsize=12,weight="bold")
    t=(f"Superficie en el parque: {d['km2']:.1f} km²\n\n"
       f"Alta:   {d['alta']:.1f} km²  ({d['alta_pct']:.1f} %)\nMedia:  {d['media']:.1f} km²  ({d['media_pct']:.1f} %)\nBaja:   {d['baja']:.1f} km²  ({d['baja_pct']:.1f} %)\n\n"
       f"Índice medio (0–1): {d['indice']:.3f}\nAltitud media: {d['elev']:.0f} m\nPendiente media: {d['pend']:.1f}°\nNDVI medio (2023–24): {d['ndvi']:.3f}")
    fig.text(.57,y0+.32,t,fontsize=10.5,va="top",family="monospace")
fig.legend(handles=[Patch(color="#bd0026",label="Alta"),Patch(color="#fd8d3c",label="Media"),Patch(color="#ffeda0",label="Baja"),Patch(fc="none",ec="#006400",ls="--",label="Límite PN (CONANP)"),Patch(fc="none",ec="#0047ab",label="Porción municipal")],loc="upper center",bbox_to_anchor=(.5,.938),ncol=5,fontsize=9,frameon=False)
foot=("Método: índice ponderado de 5 sub-índices normalizados dentro del parque (vegetación/NDVI degradado 25 %, déficit de humedad del suelo 20 %, pendiente 20 %, "
 "textura de suelo arena−arcilla 20 %, recurrencia de incendios 15 %), ~100 m. Clases alta/media/baja = terciles calculados sobre Santiago + Arteaga juntos.\n"
 "Fuentes: CONANP ANP federales (sep. 2026); geoBoundaries MEX ADM2 (World Bank/INEGI 2012); Copernicus DEM GLO-30; MODIS MOD13Q1 v061 (NDVI 2023–24); "
 "MODIS MCD64A1 v061 (área quemada 2012–2025); SoilGrids 2.0 (ISRIC); TerraClimate (humedad del suelo 2012–2021).\n"
 f"Limitaciones: según CONANP el parque está solo en Nuevo León; la porción de Arteaga (Coahuila) mide apenas ~{st['Arteaga']['km2']:.0f} km² y probablemente es un desajuste de los límites municipales de 2012 (no es el Marco Geoestadístico INEGI vigente). "
 "Los pesos son un supuesto experto sin calibrar; TerraClimate (~4 km) y MCD64A1 (500 m) son gruesos; no incluye geología/karst ni lluvia. Prototipo para priorizar; requiere validación en campo.")
import textwrap
fig.text(.05,.015,"\n".join(textwrap.fill(l,150) for l in foot.split("\n")),fontsize=6.6,va="bottom")
fig.savefig(f"{OP}/ficha_santiago_arteaga.pdf"); fig.savefig(f"{OP}/ficha_santiago_arteaga.png",dpi=130)
print(json.dumps(st,indent=1))

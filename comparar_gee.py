import sys,json,numpy as np,rasterio,geopandas as gpd
from rasterio.warp import reproject,Resampling
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap; from matplotlib.patches import Patch
P=rasterio.open("out/parque/clases_santiago_arteaga.tif"); py=P.read(1)
GIN=sys.argv[1] if len(sys.argv)>1 else "../resultados/clases_gee_v2.tif"; OUT=sys.argv[2] if len(sys.argv)>2 else "out/parque/comparacion_gee_python"; LBL=sys.argv[3] if len(sys.argv)>3 else "GEE"
G=rasterio.open(GIN); ge=np.zeros_like(py)
reproject(G.read(1),ge,src_transform=G.transform,src_crs=G.crs,dst_transform=P.transform,dst_crs=P.crs,resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
from indice_infiltracion import S,N,R
pix=(R*111.32)**2*np.cos(np.radians((S+N)/2))
gpix=(G.res[0]*111.32)**2*np.cos(np.radians(25.38)); graw=G.read(1)
out={"gee_nativo":{"crs":str(G.crs),"res_deg":G.res[0],"shape":G.shape,"nodata_declarado":G.nodata,"pix_validos":int((graw>0).sum()),"pct_validos":round(100*float((graw>0).mean()),1),
     "km2_por_clase_nativo":{n:round(float((graw==c).sum()*gpix),1) for c,n in [(3,"alta"),(2,"media"),(1,"baja")]}}}
out["km2_malla_comun"]={"gee":{n:round(float((ge==c).sum()*pix),1) for c,n in [(3,"alta"),(2,"media"),(1,"baja")]},"python":{n:round(float((py==c).sum()*pix),1) for c,n in [(3,"alta"),(2,"media"),(1,"baja")]}}
both=(py>0)&(ge>0)
out["cobertura_km2"]={"ambos":round(float(both.sum()*pix),1),"solo_python":round(float(((py>0)&(ge==0)).sum()*pix),1),"solo_gee":round(float(((ge>0)&(py==0)).sum()*pix),1)}
cm=np.array([[((py==i)&(ge==j)&both).sum() for j in (1,2,3)] for i in (1,2,3)])
n=cm.sum(); po=np.trace(cm)/n; pe=(cm.sum(1)*cm.sum(0)).sum()/n**2; k=(po-pe)/(1-pe)
adj=(np.abs(np.subtract.outer([1,2,3],[1,2,3]))<=1)
out["matriz_confusion_(filas_python_baja,media,alta;_cols_gee)"]=cm.tolist()
out["acuerdo_pct"]=round(100*po,1); out["kappa"]=round(float(k),3); out["desacuerdo_extremo_pct(alta<->baja)"]=round(100*(cm[0,2]+cm[2,0])/n,1)
d=py.astype(int)-ge.astype(int)
cat=np.zeros(py.shape,np.uint8); cat[both&(d==0)]=1; cat[both&(d==-1)]=2; cat[both&(d<=-2)]=3; cat[both&(d==1)]=4; cat[both&(d>=2)]=5; cat[(py>0)&(ge==0)]=6; cat[(ge>0)&(py==0)]=7
np.save(OUT+"_cat.npy",cat); json.dump(out,open(OUT+".json","w"),indent=1,ensure_ascii=False); print(json.dumps(out,indent=1,ensure_ascii=False))
cols=["#d9d9d9","#6baed6","#08519c","#fdae6b","#a63603","#74c476","#9e9ac8"]
labs=["Misma clase","GEE 1 clase más alta","GEE 2 clases más alta (baja→alta)","Python 1 clase más alta","Python 2 clases más alta","Solo en Python (fuera del AOI GEE)","Solo en GEE (fuera del AOI Python)"]
park=gpd.read_file("out/parque/pn_cumbres_conanp.geojson"); mun=gpd.read_file("out/parque/municipios_parque.geojson")
rows,cs=np.where(cat>0); T=P.transform
x0,x1=T.c+cs.min()*R-.02,T.c+cs.max()*R+.02; y0,y1=T.f-rows.max()*R-.02,T.f-rows.min()*R+.02
fig,ax=plt.subplots(figsize=(12,9)); ext=[T.c,T.c+P.width*R,T.f-P.height*R,T.f]
ax.imshow(np.ma.masked_equal(cat,0),cmap=ListedColormap(cols),vmin=1,vmax=7,extent=ext,interpolation="nearest")
park.boundary.plot(ax=ax,color="#006400",lw=1.5,ls="--"); mun[mun.municipio.isin(["Santiago","Arteaga"])].boundary.plot(ax=ax,color="k",lw=1)
ax.set_xlim(x0,x1); ax.set_ylim(y0,y1); ax.set_xlabel("Longitud"); ax.set_ylabel("Latitud")
ax.set_title(f"Comparación de clases {LBL} vs Python — Santiago + Arteaga (PN Cumbres)\nAcuerdo {out['acuerdo_pct']} % · kappa {out['kappa']} (en el área común, {out['cobertura_km2']['ambos']} km²)")
ax.legend(handles=[Patch(color=c,label=l) for c,l in zip(cols,labs)]+[Patch(fc="none",ec="#006400",ls="--",label="PN Cumbres (CONANP)"),Patch(fc="none",ec="k",label="Municipios 2012 (geoBoundaries)")],loc="lower left",fontsize=9)
plt.tight_layout(); plt.savefig(OUT+".png",dpi=140)

// Reto 1 · Infiltración — v2: replica el flujo de Python con los polígonos subidos como assets.
// IMPORTANTE: verifica las rutas en la pestaña "Assets" (clic en el asset > copia el "Table ID").
var PROY = 'projects/hackaton-2026-mty/assets/';
var parque = ee.FeatureCollection(PROY + 'pn_cumbres_conanp');      // CONANP shp_anp_20260924
var aoiFC  = ee.FeatureCollection(PROY + 'aoi_santiago_arteaga');   // Santiago+Arteaga ∩ parque (geoBoundaries 2012)
var munis  = ee.FeatureCollection(PROY + 'municipios_parque');      // porciones municipales dentro del parque
var PG = parque.geometry(), AOI = aoiFC.geometry();

var W = {vegetacion: 0.25, humedad: 0.20, pendiente: 0.20, suelo: 0.20, incendios: 0.15};
var SLOPE_OPT = 15, SLOPE_MAX = 40, ESCALA = 100;
var PROJ = ee.Projection('EPSG:4326').atScale(ESCALA);
// Malla idéntica a la de Python (0.001°, origen -100.8 / 26.0)
var CRS_T = [0.001, 0, -100.8, 0, -0.001, 26.0];

// --- Capas (mismas fuentes y fechas que Python) ---
var dem = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').filterBounds(PG).mosaic()
            .setDefaultProjection('EPSG:4326', null, 30)
            .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024})
            .reproject({crs: 'EPSG:4326', crsTransform: CRS_T});          // DEM promediado a ~100 m como en Python
var slope = ee.Terrain.slope(dem);
var ndvi = ee.ImageCollection('MODIS/061/MOD13Q1').filterDate('2023-01-01','2024-12-31')
            .select('NDVI').median().multiply(0.0001);
var fuego = ee.ImageCollection('MODIS/061/MCD64A1').filterDate('2012-01-01','2025-06-30')
            .select('BurnDate').map(function(i){return i.gt(0).unmask(0);}).sum();
var humedad = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate('2012-01-01','2021-12-31')
            .select('soil').mean().multiply(0.1);
var arena = ee.Image('projects/soilgrids-isric/sand_mean').select(['sand_0-5cm_mean','sand_5-15cm_mean','sand_15-30cm_mean']).reduce(ee.Reducer.mean()).divide(10);
var arcilla = ee.Image('projects/soilgrids-isric/clay_mean').select(['clay_0-5cm_mean','clay_5-15cm_mean','clay_15-30cm_mean']).reduce(ee.Reducer.mean()).divide(10);

// --- Normalización con percentiles 2–98 sobre TODO el parque, a 100 m (como Python) ---
function pcts(img, ps){
  return img.rename('v').reduceRegion({reducer: ee.Reducer.percentile(ps), geometry: PG,
    crs: 'EPSG:4326', crsTransform: CRS_T, maxPixels: 1e10});
}
function norm(img){
  var p = pcts(img, [2, 98]);
  var lo = ee.Number(p.get('v_p2')), hi = ee.Number(p.get('v_p98'));
  return img.subtract(lo).divide(hi.subtract(lo)).clamp(0, 1);
}
var sVeg = ee.Image(1).subtract(norm(ndvi));
var sHum = ee.Image(1).subtract(norm(humedad));
var sPen = ee.Image(1).subtract(slope.subtract(SLOPE_OPT).divide(SLOPE_MAX - SLOPE_OPT)).clamp(0, 1);
var fmax = ee.Number(fuego.rename('v').reduceRegion({reducer: ee.Reducer.max(), geometry: PG,
             crs: 'EPSG:4326', crsTransform: CRS_T, maxPixels: 1e10}).get('v')).max(1);
var sFue = fuego.divide(fmax).clamp(0, 1);                       // meses quemados / máximo del parque
var sSue = norm(arena.subtract(arcilla));

var idx = sVeg.multiply(W.vegetacion).add(sHum.multiply(W.humedad)).add(sPen.multiply(W.pendiente))
          .add(sSue.multiply(W.suelo)).add(sFue.multiply(W.incendios))
          .divide(W.vegetacion + W.humedad + W.pendiente + W.suelo + W.incendios)
          .reproject({crs: 'EPSG:4326', crsTransform: CRS_T})
          .clip(PG).rename('indice');

// --- Terciles SOLO sobre el AOI Santiago+Arteaga (como Python) ---
var t = idx.reduceRegion({reducer: ee.Reducer.percentile([33.333, 66.667]), geometry: AOI,
          crs: 'EPSG:4326', crsTransform: CRS_T, maxPixels: 1e10});
var t1 = ee.Number(t.get('indice_p33')), t2 = ee.Number(t.get('indice_p67'));
var clase = ee.Image(1).where(idx.gte(t1), 2).where(idx.gte(t2), 3)
             .updateMask(idx.mask()).clip(AOI).rename('clase').toByte();   // 1 baja, 2 media, 3 alta

Map.centerObject(AOI, 11);
Map.addLayer(idx.clip(AOI), {min: 0.2, max: 0.8, palette: ['ffffcc','fd8d3c','bd0026']}, 'Índice continuo', false);
Map.addLayer(clase, {min: 1, max: 3, palette: ['ffeda0','fd8d3c','bd0026']}, 'Prioridad: baja/media/alta');
Map.addLayer(ee.Image().paint(parque, 0, 2), {palette: ['006400']}, 'PN Cumbres (CONANP)');
Map.addLayer(ee.Image().paint(munis, 0, 1), {palette: ['444444']}, 'Municipios (porción en parque)');
Map.addLayer(ee.Image().paint(aoiFC, 0, 3), {palette: ['0047ab']}, 'AOI Santiago + Arteaga');
print('Umbrales terciles (Python: 0.3653 / 0.5080)', t);
var areas = ee.Image.pixelArea().divide(1e6).addBands(clase).reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'clase'}), geometry: AOI,
  crs: 'EPSG:4326', crsTransform: CRS_T, maxPixels: 1e10});
print('Área km² por clase (1 baja, 2 media, 3 alta)', areas);

Export.image.toDrive({image: clase.unmask(0).toByte(), description: 'clases_gee_v2', fileNamePrefix: 'clases_gee_v2',
  folder: 'EarthEngine', region: AOI.bounds(), crs: 'EPSG:4326', crsTransform: CRS_T, maxPixels: 1e10});

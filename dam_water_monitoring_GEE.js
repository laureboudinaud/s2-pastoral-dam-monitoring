// =============================================================================
// Title      : Water Presence and Seasonal Persistence Metrics derived from
//              Sentinel-2 Time Series for Dam Monitoring in Semi-Arid Regions
// Author     : Laure Boudinaud
// Date       : 2026
// License    : CC BY 4.0
// -----------------------------------------------------------------------------
// Description:
//   Computes annual Sentinel-2 composites (2016–2025) per dam buffer zone and
//   exports zonal statistics as three long-format CSVs (one row per dam per year)
//   for post-processing in Python.
//
//   COMPOSITES (per year):
//     ndwi_p90   Annual 90th-percentile NDWI        — inter-annual water presence
//     wet_p90    Wet season (May–Oct) p90 NDWI      — peak water signal
//     dry_p90    Dry season (Jan–Feb) p90 NDWI      — residual dry-season water
//     spi        wet_p90 − dry_p90, water-masked    — Seasonal Persistence Index
//                low SPI  = year-round retention (permanent)
//                high SPI = dry-season loss (seasonal/ephemeral)
//     water      Binary mask: ndwi_p90 > 0          — water presence flag
//
//   EXPORTED TABLES:
//     dam_stats_mean_2016_2025.csv     — mean ndwi_p90, wet_p90, dry_p90, spi
//                                         water_fraction (mean of binary mask)
//     dam_stats_max_2016_2025.csv      — max pixel values (less land dilution)
//     dam_stats_spi_stddev_2016_2025.csv — SPI stddev within water surface
//                                          high stddev = permanent core + seasonal
//                                          fringe → restoration-relevant
// =============================================================================


// =============================================================================
// 1. USER PARAMETERS
// =============================================================================
var dams_pts = ee.FeatureCollection('projects/ee-laureboudinaud/assets/dam_locations_civ_tchbou');

var BUFFER_RADIUS = 1000;         // metres

var YEAR_START = 2016;
var YEAR_END   = 2025;

// Wet season (northern Côte d'Ivoire: May–October)
var WET_START = 5;
var WET_END   = 10;

// Peak dry season (January–February)
var DRY_START = 1;
var DRY_END   = 2;

// Snapshot years for map display and raster export
var SNAPSHOT_YEARS = [2016, 2021, 2025];

// Visualisation params
var NDWI_VIZ = {min: -0.4, max: 0.3, palette: ['brown', 'white', 'blue']};
var SPI_VIZ  = {min: 0,     max: 0.5, palette: ['blue', 'white', 'yellow']};

// Google Drive export folder
var EXPORT_FOLDER = 'GEE';


// =============================================================================
// 2. SETUP
// =============================================================================
var years       = ee.List.sequence(YEAR_START, YEAR_END);
var dams_buffer = dams_pts.map(function(f) { return f.buffer(BUFFER_RADIUS); });
var aoi         = dams_buffer.geometry().bounds();

Map.centerObject(dams_buffer, 9);
Map.addLayer(
  dams_buffer.style({color: 'ffffff', fillColor: '00000000', width: 1.5}),
  {}, 'Dam buffers', true
);

function maskS2clouds(img) {
  var qa            = img.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return img.updateMask(mask).copyProperties(img, ['system:time_start']);
}

function addNDWI(img) {
  return img
    .addBands(img.normalizedDifference(['B3', 'B8']).rename('NDWI'))
    .copyProperties(img, ['system:time_start']);
}

function s2(year) {
  year = ee.Number(year);
  return ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
    .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year, 12, 31))
    .filterBounds(dams_buffer)
    .map(maskS2clouds)
    .map(addNDWI);
}


// =============================================================================
// 3. ANNUAL COMPOSITES
// =============================================================================
// For each year, compute:
//   - annual_p90  : annual p90 NDWI (all months)
//   - wet_p90     : wet season p90 NDWI
//   - dry_p90     : dry season p90 NDWI
//   - spi         : wet_p90 - dry_p90

var annualMetrics = years.map(function(year) {
  year = ee.Number(year);
  var col = s2(year);

  var ndwi_p90 = col
    .select('NDWI')
    .reduce(ee.Reducer.percentile([90]))
    .rename('ndwi_p90');

  var wet_p90 = col
    .filter(ee.Filter.calendarRange(WET_START, WET_END, 'month'))
    .select('NDWI')
    .reduce(ee.Reducer.percentile([90]))
    .rename('wet_p90');

  var dry_p90 = col
    .filter(ee.Filter.calendarRange(DRY_START, DRY_END, 'month'))
    .select('NDWI')
    .reduce(ee.Reducer.percentile([90]))
    .rename('dry_p90');

  // Binary water mask from annual p90 (for SPI masking and visualisation)
  var water = ndwi_p90.gt(0).rename('water');
  var spi = wet_p90.subtract(dry_p90).updateMask(water).rename('spi');

  return ndwi_p90
    .addBands(wet_p90)
    .addBands(dry_p90)
    .addBands(spi)
    .addBands(water)
    .set('year', year);
});

var metricsCol = ee.ImageCollection(annualMetrics);


// =============================================================================
// 4. ZONAL STATISTICS — THREE SEPARATE TABLES
// =============================================================================
// One row per dam per year.
// Columns: system:index, year, ndwi_p90, wet_p90, dry_p90, spi,
//          water_fraction (share of pixels with NDWI > 0 in the buffer)

function reduceByYear(reducer, bands) {
  return years.map(function(year) {
    year = ee.Number(year);
    var img = metricsCol.filter(ee.Filter.eq('year', year)).first();
    var stats = img
      .select(bands)
      .reduceRegions({
        collection: dams_buffer,
        reducer:    reducer,
        scale:      10,
        crs:        'EPSG:4326'
      });
    return stats.map(function(f) {
      var centroid = f.geometry().centroid(1);
      return f.set('year', year)
               .set('lon', centroid.coordinates().get(0))
               .set('lat', centroid.coordinates().get(1));
    });
  });
}

// --- Table A: MEAN (mean over buffer) ---
// ndwi_p90, wet_p90, dry_p90, spi, water (mean of binary = water_fraction)
var meanStats = ee.FeatureCollection(
  reduceByYear(
    ee.Reducer.mean(),
    ['ndwi_p90', 'wet_p90', 'dry_p90', 'spi', 'water'],
    'mean'
  )
).flatten();

// --- Table B: MAX (peak pixel value within buffer) ---
// ndwi_p90, wet_p90, dry_p90, spi — 
// (captures the dam signal, less diluted by surrounding land)
// Note: max SPI useful for identifying the most seasonal pixel in the buffer
var maxStats = ee.FeatureCollection(
  reduceByYear(
    ee.Reducer.max(),
    ['ndwi_p90', 'wet_p90', 'dry_p90', 'spi'],
    'max'
  )
).flatten();

// --- Table C: STDDEV of SPI — computed over water-detected pixels only ---
// high stddev = permanent core + seasonal fringe within the water surface
//             = spatially heterogeneous retention → restoration-relevant
var stddevStats = ee.FeatureCollection(
  reduceByYear(
    ee.Reducer.stdDev().setOutputs(['spi_stddev']),
    ['spi'],
    'stddev'
  )
).flatten();


// =============================================================================
// 5. MAP — SNAPSHOT YEARS
// =============================================================================
// Shows p90 NDWI and SPI (masked to water pixels) for three years.
// Allows visual inspection of change before further analysis.

SNAPSHOT_YEARS.forEach(function(yr) {
  var col = s2(yr);

  var ndwi_p90 = col.select('NDWI')
    .reduce(ee.Reducer.percentile([90])).rename('ndwi_p90');

  var wet_p90 = col
    .filter(ee.Filter.calendarRange(WET_START, WET_END, 'month'))
    .select('NDWI').reduce(ee.Reducer.percentile([90]));

  var dry_p90 = col
    .filter(ee.Filter.calendarRange(DRY_START, DRY_END, 'month'))
    .select('NDWI').reduce(ee.Reducer.percentile([90]));

  var water = ndwi_p90.gt(0).selfMask();
  var spi   = wet_p90.subtract(dry_p90).updateMask(water).rename('spi');

  Map.addLayer(ndwi_p90.clip(dams_buffer), NDWI_VIZ, 'NDWI p90 ' + yr, false);
  Map.addLayer(spi.clip(dams_buffer), SPI_VIZ, 'SPI ' + yr, false);
  Map.addLayer(water.clip(dams_buffer), {palette: 'blue'}, 'Water mask ' + yr, false);
});


// =============================================================================
// 6. EXPORT
// =============================================================================
// Primary outputs for post-processing and further analysis in Python or R
Export.table.toDrive({
  collection:  meanStats,
  description: 'dam_stats_mean_2016_2025',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

Export.table.toDrive({
  collection:  maxStats,
  description: 'dam_stats_max_2016_2025',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

Export.table.toDrive({
  collection:  stddevStats,
  description: 'dam_stats_spi_stddev_2016_2025',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

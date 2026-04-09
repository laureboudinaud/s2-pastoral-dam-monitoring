// =============================================================================
// Title      : Sentinel-2 NDWI Time Series and Seasonal Persistence Index (SPI)
//              for Pastoral Dam Functionality and Water Availability Assessment
// Author     : Laure Boudinaud
// Date       : February 2026
// Description: Computes annual NDWI percentile composites and a Seasonal
//              Persistence Index (SPI) from Sentinel-2 imagery to assess
//              two dimensions of dam water availability for pastoral
//              infrastructure monitoring:
//                (1) inter-annual functionality — is water present across years?
//                (2) intra-annual persistence  — does water persist through
//                    the dry season within a given year?
//              Developed for the Bounkani and Tchologo regions of northern
//              Côte d'Ivoire; adaptable to any region with a defined area
//              of interest and dam point layer.
// License    : CC BY 4.0
// =============================================================================


// =============================================================================
// 1. USER PARAMETERS
// =============================================================================
// -- Area of interest --
// Note: civ_adm1 is a private asset. Replace with your own boundary layer.
// See README for alternative public boundary options.
var aoi = civ_adm1.filter(ee.Filter.inList("ADM1_FR", ["Bounkani", "Tchologo"]));

// -- Dam point layer --
var dams_pts = ee.FeatureCollection("users/your_username/your_dam_asset")
                 .filterBounds(aoi);

// -- Buffer radius around dam points (in metres) --
// 1000 m recommended for pastoral dams in northen Côte d'Ivoire
var BUFFER_RADIUS = 1000;

// -- Reference years for interannual comparison --
var YEAR_REF  = 2018;   // baseline  year
var YEAR_CURR = 2025;   // analysis year

// -- Wet season months (inclusive) --
// Adjust to match the rainy season of your study area (northern Côte d'Ivoire: May–October)
var WET_START = 5;   // May
var WET_END   = 10;  // October

// -- Dry season months (inclusive) --
// Adjust to match the peak dry season of your study area (northern Côte d'Ivoire: January–February)
var DRY_START = 1;   // January
var DRY_END   = 2;   // February

// -- SPI visualization range --
// SPI = wet_NDWI_p90 - dry_NDWI_p90
var SPI_VIZ_MIN = 0; // Low SPI  : permanent/stable water (present in both seasons)
var SPI_VIZ_MAX = 0.5; // High SPI : strongly seasonal water (present in wet, absent in dry)

// -- Output folder on Google Drive for exports --
var EXPORT_FOLDER = 'GEE';


// =============================================================================
// 2. SETUP
// =============================================================================
Map.centerObject(aoi, 9);

// Create buffer zones around each dam point
var bufferFn = function(f) { return f.buffer(BUFFER_RADIUS); };
var dams_buffer = dams_pts.map(bufferFn);


// =============================================================================
// 3. FUNCTIONS
// =============================================================================
// -- Cloud masking using Sentinel-2 QA60 band --
function maskS2clouds(img) {
  var qa = img.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return img.updateMask(mask)
            .copyProperties(img, ['system:time_start']);
}

// -- Load Sentinel-2 collection and compute NDWI for a given year --
var s2NDWI = function(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = ee.Date.fromYMD(year, 12, 31);
  return ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
            .filterDate(start, end)
            .filterBounds(aoi)
            .map(maskS2clouds)
            .map(function(image) {
              var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
              return image.addBands(ndwi)
                          .copyProperties(image, ['system:time_start']);
            });
};

// -- Compute 90th percentile NDWI composite --
var p90 = function(collection) {
  return collection.select('NDWI')
                   .reduce(ee.Reducer.percentile([90]))
                   .rename('ndwi_p90');
};

// -- Season filter functions --
var wetSeason = function(collection) {
  return collection.filter(ee.Filter.calendarRange(WET_START, WET_END, 'month'));
};
var drySeason = function(collection) {
  return collection.filter(ee.Filter.calendarRange(DRY_START, DRY_END, 'month'));
};


// =============================================================================
// 4. COMPUTE NDWI AND WATER EXTENT
// =============================================================================
var all_ref  = s2NDWI(YEAR_REF);
var all_curr = s2NDWI(YEAR_CURR);

// Annual p90 NDWI composites (full year, all seasons combined)
var ndwi_p90_ref  = p90(all_ref);
var ndwi_p90_curr = p90(all_curr);

// Binary water masks: pixels with NDWI > 0 are classified as water
var water_ref  = ndwi_p90_ref.gt(0).selfMask();
var water_curr = ndwi_p90_curr.gt(0).selfMask();


// =============================================================================
// 5. COMPUTE SEASONAL PERSISTENCE INDEX (SPI)
// =============================================================================
// Seasonal p90 NDWI composites clipped to dam buffer zones
var wet_ref_p90  = p90(wetSeason(all_ref)).clip(dams_buffer);
var dry_ref_p90  = p90(drySeason(all_ref)).clip(dams_buffer);
var wet_curr_p90 = p90(wetSeason(all_curr)).clip(dams_buffer);
var dry_curr_p90 = p90(drySeason(all_curr)).clip(dams_buffer);

// SPI = NDWI_p90_wet - NDWI_p90_dry
var spi_ref  = wet_ref_p90.subtract(dry_ref_p90).rename('SPI_' + YEAR_REF);
var spi_curr = wet_curr_p90.subtract(dry_curr_p90).rename('SPI_' + YEAR_CURR);

// Mask SPI to pixels where water was detected in the current year.
var spi_curr_masked = spi_curr.updateMask(water_curr);


// =============================================================================
// 6. VISUALIZATION
// =============================================================================
var ndwi_viz = {min: -0.43, max: 0.1, palette: ['brown', 'white', 'blue']};
var spi_viz  = {min: SPI_VIZ_MIN, max: SPI_VIZ_MAX, palette: ['blue', 'white', 'yellow']};

// Annual NDWI p90 composites
Map.addLayer(ndwi_p90_ref.clip(dams_buffer),
             ndwi_viz, 'NDWI p90 ' + YEAR_REF + ' (dam buffers)');
Map.addLayer(ndwi_p90_curr.clip(dams_buffer),
             ndwi_viz, 'NDWI p90 ' + YEAR_CURR + ' (dam buffers)');

// Binary water masks
Map.addLayer(water_ref.clip(dams_buffer),
             {palette: '#ADD8E6'}, 'Water presence ' + YEAR_REF, false);
Map.addLayer(water_curr.clip(dams_buffer),
             {palette: '0000FF'}, 'Water presence ' + YEAR_CURR, false);

// Seasonal Persistence Index
Map.addLayer(spi_curr_masked.clip(dams_buffer),
             spi_viz, 'SPI ' + YEAR_CURR +
             ' (blue=permanent, yellow=seasonal)');

// Dam buffer outlines for reference
Map.addLayer(dams_buffer.style({color: 'ffffff', fillColor: '00000000', width: 1.5}),
             {}, 'Dam buffers', false);


// =============================================================================
// 7. EXPORT
// =============================================================================
// -- Export mean SPI value per dam buffer as CSV --
var spi_by_buffer = spi_curr_masked.reduceRegions({
  collection: dams_buffer,
  reducer:    ee.Reducer.mean(),   // pixels outside water mask are excluded
  scale:      10,
  crs:        'EPSG:4326'
});

Export.table.toDrive({
  collection:  spi_by_buffer,
  description: 'SPI_' + YEAR_CURR + '_by_dam_buffer',
  folder:      EXPORT_FOLDER,
  fileFormat:  'CSV'
});

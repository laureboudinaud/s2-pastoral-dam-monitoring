# s2-pastoral-dam-monitoring
GEE script and datasets for assessing dam water availability in semi-arid regions (northern Côte d'Ivoire as example study area) using Sentinel-2 NDWI time series and a Seasonal Persistence Index (SPI)

**Author:** Laure Boudinaud  
**Date:** February 2026  
**License:** CC BY 4.0  
**Associated publication:** Boudinaud, L. et al. (2026) *(in preparation)*

\---

## Overview

This repository contains a Google Earth Engine (GEE) JavaScript script to assess two dimensions of pastoral dam water availability using Sentinel-2 time series:

1. **Inter-annual functionality** — whether surface water is detectable at a dam location across years
2. **Intra-annual persistence** — whether water is maintained through the dry season within a given year

The script was developed for the Bounkani and Tchologo regions of northern Côte d'Ivoire and is fully adaptable to any region with a defined area of interest and dam point layer.

\---

## Repository Contents

|File|Description|
|-|-|
|`dam\_water\_monitoring\_GEE.js`|Main GEE script|
|`README.md`|This file|
|`SPI\_2025\_by\_dam\_buffer.csv`|Exported SPI values per dam buffer|
|`dam\_locations\_EO\_classified.csv`|Dam locations with EO-derived functionality classifications|

\---

## Requirements

* A Google Earth Engine account : https://earthengine.google.com
* A dam point layer uploaded as a GEE asset (see Input Data below)
* A boundary layer for your area of interest (see Area of Interest below)

\---

## How to Run

1. Open the GEE Code Editor at https://code.earthengine.google.com
2. Create a new script and paste the contents of `dam\_water\_monitoring\_GEE.js`
3. Update the USER PARAMETERS in Section 1 of the script
4. Click **Run**

\---

## Input Data

### Dam Point Layer

The script requires a point FeatureCollection representing dam locations. Each point should correspond to the centroid or spillway location of a dam.



**Study area dataset:** For the Bounkani and Tchologo regions, the dam point layer was built in two stages. An initial set of 68 points was derived from field surveys and validated against Sentinel-2 imagery. This set was then expanded using the water detection component of this script (Section 4) to identify additional water bodies consistent with dam infrastructure, reaching a total of 312 dam locations.



**To use your own dam layer:**

1. Prepare a shapefile or CSV with point coordinates
2. Upload it to GEE via: **Assets > New > Shapefiles or CSV**
3. Replace the asset path in the USER PARAMETERS section of the script:

```javascript
var dams\\\_pts = ee.FeatureCollection("users/your\\\_username/your\\\_dam\\\_asset")
                 .filterBounds(aoi);
```

OR

**Deriving dam points directly from the script:**  
If no field-surveyed point layer is available, candidate dam locations can be extracted directly from the binary water mask computed in Section 4 of the script. To do this, convert the `water\_curr` image to vectors and extract centroids:

```javascript
// Extract candidate dam points from detected water polygons
var water\_polygons = water\_curr.reduceToVectors({
  geometry: aoi,
  scale: 10,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'water',
  crs: 'EPSG:4326'
});

var candidate\_points = water\_polygons.map(function(f) {
  return f.centroid();
});
```

> Important: Automatically derived points should be manually reviewed before use. The water detection step identifies all surface water, including rivers, flooded areas, and other non-dam features. Visual inspection against high-resolution imagery (e.g. Google Satellite in the GEE basemap) is strongly recommended before using these points for dam functionality assessment.



### Area of Interest

The script as provided references a private administrative boundary asset (`civ\_adm1`) specific to Côte d'Ivoire. To adapt to a different study area, replace this with one of the following publicly available options in GEE:

|Dataset|Asset path|Notes|
|-|-|-|
|FAO GAUL 2015|`FAO/GAUL/2015/level1`|Good global coverage but predates recent administrative reorganisations in some countries|
|Custom upload|`users/your\_username/your\_boundary\_asset`|Recommended when public datasets do not reflect recent administrative subdivisions — sources may be FAO GAUL (2024), GeoBoundaries, OCHA/HDX|

> \*\*Note for Côte d'Ivoire users:\*\* The Bounkani and Tchologo regions are recent administrative subdivisions that postdate the FAO GAUL 2015 dataset. Users working in this area should upload a current boundary layer from \[OCHA HDX](https://data.humdata.org) or the Ivorian national mapping agency (BNETD-CIGN).

\---

## User Parameters

All parameters are consolidated in **Section 1** of the script:

|Parameter|Description|Default|
|-|-|-|
|`aoi`|Area of interest boundary|Bounkani \& Tchologo, Côte d'Ivoire|
|`dams\_pts`|Dam point layer|Private asset — replace with your own|
|`BUFFER\_RADIUS`|Buffer radius around dam points (metres)|`1000`|
|`YEAR\_REF`|Baseline / reference year|`2018`|
|`YEAR\_CURR`|Current / analysis year|`2025`|
|`WET\_START`|Start month of wet season (inclusive)|`5` (May)|
|`WET\_END`|End month of wet season (inclusive)|`10` (October)|
|`DRY\_START`|Start month of peak dry season (inclusive)|`1` (January)|
|`DRY\_END`|End month of peak dry season (inclusive)|`2` (February)|
|`SPI\_VIZ\_MIN`|Minimum value for SPI visualisation|`0`|
|`SPI\_VIZ\_MAX`|Maximum value for SPI visualisation|`0.5`|
|`EXPORT\_FOLDER`|Google Drive folder name for exported outputs|`GEE`|

\---

## Outputs

### Map Layers (visualised in the GEE Code Editor)

* Annual NDWI p90 composites for reference and current year, clipped to dam buffers
* Binary water presence masks for both years
* Seasonal Persistence Index (SPI) for the current year, masked to water-detected pixels
* Dam buffer outlines

### Exported Files (saved to Google Drive)

* `SPI\_\[YEAR]\_by\_dam\_buffer.csv` — Mean SPI value per dam buffer, used to classify dams into water persistence categories

### SPI Classification

|Mean SPI|Interpretation|
|-|-|
|< 0.1|Permanent water|
|0.1 – 0.3|Semi-permanent|
|> 0.3|Seasonal / ephemeral|
|No value|No water detected (likely non-functional dam)|

\---

## Citation

If you use this script in your research, please cite both the associated publication and the code deposit:

**Publication:**

> Boudinaud, L. et al. (2026) \*(in preparation)\*

**Code:**

> Boudinaud, L. (2026): GEE script for Sentinel-2 NDWI and Seasonal Persistence Index for pastoral dam monitoring, Côte d'Ivoire, \*(in preparation)\*

\---

## Contact

Laure Boudinaud  
For questions or issues, please contact laure@alateo-analytics.com or laure.boudinaud@gmail.com

\---

## Acknowledgements

This script was developed in the framework of a project on agro-pastoral resource monitoring and conflict dynamics in northern Côte d'Ivoire, implemented by PHI Consulting in partnership with the Réseau Billital Maroobé (RBM) and the International Organization for Migration (IOM), with funding from the Government of Japan.


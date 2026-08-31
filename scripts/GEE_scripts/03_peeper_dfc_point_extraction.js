// =======================================================================================
// Extract predicted Date of First Calling (DFC) at NAAMP route centroids.
// Values are computed at points only, so this is not limited by the New York
// extent used for the gridded app asset.
//
// NOTE: run in 8-9 year chunks. Change chunkStartYear / chunkEndYear per run.
// =======================================================================================
 
// ===== ROUTE POINTS =====
var rawRoutes = ee.FeatureCollection('projects/predicted-peeper-calling-onset/assets/naamp_route_coords');
 
// CSV upload gives a table without geometry, so build points from the columns
var routes = rawRoutes.map(function(feature) {
  var lon = feature.get('lon');
  var lat = feature.get('lat');
  return feature.setGeometry(ee.Geometry.Point([lon, lat]));
});
 
// ===== DATA LOAD =====
var era5LandHourly = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY');
 
// ===== DEFINE CONSTANTS FOR CALCULATION =====
// TS3 threshold from Lovett (2013)
var ts3Threshold = 44.3;
var baseTemp = 3;
var searchStartMonth = 2; // February 1
var searchStartDay = 1;
var searchDays = 150; // search through end of June
 
// ===== CONFIGURE CHUNKS =====
// GEE can only handle a handful of years at a time with hourly ERA5.
// Start with 9 years; narrow further if tasks fail.
var chunkStartYear = 2006;
var chunkEndYear = 2015;
 
// ===== LOGIC TO GET PREDICTED DFC =====
// Following Lovett (2013): accumulate GDD above 3°C from Feb 1
// Return the day of year when cumulative TS3 first exceeds 44.3
var getPredictedDFC = function(year) {
  year = ee.Number(year);
 
  var startDate = ee.Date.fromYMD(year, searchStartMonth, searchStartDay);
  var days = ee.List.sequence(0, searchDays - 1);
 
  // filter the whole window once, up front
  var hourly = era5LandHourly
    .filterDate(startDate, startDate.advance(searchDays, 'day'))
    .select('temperature_2m');
 
  var dailyGDD = ee.ImageCollection(days.map(function(dayOffset) {
    var d = startDate.advance(ee.Number(dayOffset), 'day');
    return hourly.filterDate(d, d.advance(1, 'day'))
      .mean()
      .subtract(273.15) // convert to Celsius
      .subtract(baseTemp)
      .max(0) // any temperature below 0 is considered 0
      .set('system:time_start', d.millis());
  }));
 
  // create array of cumulative sum for each day in each pixel
  var cumSum = dailyGDD.sort('system:time_start').toArray().arrayAccum(0);
 
  // count days still below threshold '1' = index of first crossing
  var nBelow = cumSum.lt(ts3Threshold)
    .arrayReduce(ee.Reducer.sum(), [0])
    .arrayProject([0])
    .arrayFlatten([['n']]);
 
  var predictedDFC = nBelow.add(32)
    .rename('predicted_dfc')
    .set({'year': year, 'system:time_start': startDate.millis()});
 
  return predictedDFC.updateMask(nBelow.lt(searchDays));
};
 
// ===== EXTRACT AT ROUTE POINTS, ONE ROW PER ROUTE x YEAR =====
var extractYear = function(year) {
  year = ee.Number(year);
 
  return getPredictedDFC(year)
    .reduceRegions({
      collection: routes,
      reducer: ee.Reducer.first(),
      scale: 11132,
      crs: 'EPSG:4326'
    }).map(function(f) {
      return f.set('year', year);
    });
};
 
var chunkExtracted = ee.FeatureCollection(
  ee.List.sequence(chunkStartYear, chunkEndYear).map(function(year) {
    return extractYear(year);
  })
).flatten();
 
Export.table.toDrive({
  collection: chunkExtracted,
  description: 'naamp_route_predicted_dfc_' + chunkStartYear + '_' + chunkEndYear,
  fileFormat: 'CSV',
  selectors: ['RouteNumber', 'year', 'first']
});
 
// ===== SANITY CHECK BEFORE RUNNING THE EXPORT =====
// Confirm the points landed where you expect and that values come back.
print('Route count:', routes.size());
print('First route:', routes.first());
print('Sample extraction (1 year):',
  extractYear(chunkStartYear).limit(5));
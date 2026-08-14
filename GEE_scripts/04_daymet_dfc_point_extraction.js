// =======================================================================================
// Extract predicted Date of First Calling (DFC) at NAAMP route centroids
// using DAYMET V4 (1 km) instead of ERA5-Land (~9 km).
//
// Purpose: test whether ERA5's coarse resolution smooths climate variation
// that matters for predicting calling onset. Same model, same threshold,
// same sites and years - only the climate input changes.
//
// NOTE: Daymet begins in 1980, so 1997-2015 is fully covered.
// Run in chunks. Change chunkStartYear / chunkEndYear per run.
// =======================================================================================
 
// ===== ROUTE POINTS =====
var rawRoutes = ee.FeatureCollection('projects/predicted-peeper-calling-onset/assets/naamp_route_coords');
 
var routes = rawRoutes.map(function(feature) {
  var lon = feature.get('lon');
  var lat = feature.get('lat');
  return feature.setGeometry(ee.Geometry.Point([lon, lat]));
});
 
// ===== DATA LOAD =====
// Daymet supplies daily tmax and tmin (already in Celsius), not hourly temperature.
// Daily mean is therefore (tmax + tmin) / 2 rather than a 24-hour average.
var daymet = ee.ImageCollection('NASA/ORNL/DAYMET_V4');
 
// ===== DEFINE CONSTANTS FOR CALCULATION =====
// TS3 threshold from Lovett (2013)
var ts3Threshold = 44.3;
var baseTemp = 3;
var searchStartMonth = 2; // February 1
var searchStartDay = 1;
var searchDays = 150; // search through end of June
var scale = 1000;     // Daymet native resolution
 
// ===== CONFIGURE CHUNKS =====
var chunkStartYear = 2006;
var chunkEndYear = 2015;
 
// ===== LOGIC TO GET PREDICTED DFC =====
// Following Lovett (2013): accumulate GDD above 3 C from Feb 1
// Return the day of year when cumulative TS3 first exceeds 44.3
var getPredictedDFC = function(year) {
  year = ee.Number(year);
 
  var startDate = ee.Date.fromYMD(year, searchStartMonth, searchStartDay);
 
  var dailyGDD = daymet
    .filterDate(startDate, startDate.advance(searchDays, 'day'))
    .map(function(img) {
      var meanTemp = img.select('tmax').add(img.select('tmin')).divide(2);
      return meanTemp
        .subtract(baseTemp)
        .max(0) // any value below 0 contributes 0
        .copyProperties(img, ['system:time_start']);
    });
 
  // cumulative sum across the accumulation window, per pixel
  var cumSum = ee.ImageCollection(dailyGDD)
    .sort('system:time_start')
    .toArray()
    .arrayAccum(0);
 
  // count days still below threshold = offset of first crossing
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
      scale: scale,
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
  description: 'naamp_route_predicted_dfc_daymet_' + chunkStartYear + '_' + chunkEndYear,
  fileFormat: 'CSV',
  selectors: ['RouteNumber', 'year', 'first']
});
 
// ===== SANITY CHECK BEFORE RUNNING THE EXPORT =====
// Confirm values come back and are in a plausible range before queuing tasks.
var sample = extractYear(chunkStartYear);
print('Route count:', routes.size());
print('Non-null returns:', sample.filter(ee.Filter.notNull(['first'])).size());
print('Sample values:', sample.limit(10).aggregate_array('first'));
// ============================================================
// Script to Export Predicted Date of First Call (DFC) for each
// pixel in each year.
// Run once to export asset for application 
// ============================================================

// ===== EXPORT FOLDER =====
var assetRoot = 'projects/predicted-peeper-calling-onset/assets/dfc_annual';

// ===== DATA LOAD ===== 
var era5LandHourly = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY');
var peeperRange = ee.FeatureCollection('projects/predicted-peeper-calling-onset/assets/spring_peeper_range_iucn');
var states = ee.FeatureCollection('TIGER/2018/States');
var ny = states.filter(ee.Filter.eq('NAME', 'New York'));
var geom = peeperRange.geometry().intersection(ny.geometry(), 1000);


// ===== DEFINE CONSTANTS FOR CALCULATION ===== 
// year range
var currentYear = new Date().getFullYear() - 1;

// TS3 threshold from Lovett (2013)
var ts3Threshold = 44.3;
var baseTemp = 3;
var searchStartMonth = 2; // February 1
var searchStartDay = 1;
var searchDays = 150; // search through end of June

// ===== LOGIC TO GET PREDICTED DFC =====
// Function to get predicted DFC for one year
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
    .toFloat()
    .rename('predicted_dfc')
    .set({'year': year, 'system:time_start': startDate.millis()});

  return predictedDFC.updateMask(nBelow.lt(searchDays));
};

// ===== EXPORT PREDICTED DFC - ONE IMAGE PER YEAR =====
for (var y = 1950; y <= currentYear; y++) {
  Export.image.toAsset({
    image: getPredictedDFC(y),
    description: 'dfc_' + y,
    assetId: assetRoot + '/dfc_' + y,
    region: geom,
    scale: 11132,
    crs: 'EPSG:4326',
    maxPixels: 1e9
  });
}
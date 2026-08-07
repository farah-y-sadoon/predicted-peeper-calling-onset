// ============================================================
// Script to run application
// ============================================================

// ===== DATA LOAD ===== 
// Load image collection for precomputed date of first calling (DFC) 
var annualDFCCollection = ee.ImageCollection('projects/predicted-peeper-calling-onset/assets/dfc_annual');

// Load shape files for geometry
var peeperRange = ee.FeatureCollection('projects/predicted-peeper-calling-onset/assets/spring_peeper_range_iucn');
var peeperRangeGeometry = peeperRange.geometry();
var states = ee.FeatureCollection('TIGER/2018/States');
var ny = states.filter(ee.Filter.eq('NAME', 'New York'));
var geom = peeperRange.geometry().intersection(ny.geometry(), 1000);

// ===== GET BASELINE MEAN DFC (1950-1979) =====
var baselineStart = 1950;
var baselineEnd = 1979;
 
var baselineDFC = annualDFCCollection
  .filter(ee.Filter.and(
    ee.Filter.gte('year', baselineStart), 
    ee.Filter.lte('year', baselineEnd)
    ))
  .mean()
  .rename('baseline_dfc');
  
// ===== DEFINE ANOMALY COLLECTION ====   
// Create image collection with each year's deviation from the baseline (1950 - 2025)
var anomalyCollection = annualDFCCollection.map(function(image) {
  return ee.Image(image).select('predicted_dfc')
    .subtract(baselineDFC)
    .rename('dfc_anomaly')
    .copyProperties(image, ['year', 'system:time_start']);
});
 
// ===== UI FUNCTIONALITY =====
// Constants
var anomalyVis = {
  min: -35,
  max: 35,
  palette: ['#2166ac', '#67a9cf', '#f7f7f7', '#ef8a62', '#b2182b']
};
 
var layerAnomaly = 0;
var layerMarker = 1;
var anomalyLimit = 35
 
var anomalyList = anomalyCollection.sort('year').toList(76);
 
// Layout Scaffolding
ui.root.clear();
 
var map = ui.Map();
var mainPanel = ui.Panel({style: {width: '400px', padding: '12px'}});
 
// Widgets
var title = ui.Label({
  value: 'Predicted Spring Peeper (Pseudacris crucifer) Calling Onset', 
  style: {fontSize: '18px', fontWeight: 'bold', margin: '0 0 4px 0'}
});
 
var subtitle = ui.Label({
  value: 'New York State, 1950 - 2025', 
  style: {fontSize: '14px', color: '#666', margin: '0 0 12px 0'}
});
 
var methodText = ui.Label(
  'Spring peepers (Pseudacris crucifer) begin calling once accumulated ' +
  'warmth reaches a threshold. Following Lovett (2013), this tool predicts ' +
  'the date of first calling (DFC) as the day on which degree-days above 3 °C, ' +
  'accumulated from February 1, first reach 44.3.',
  {fontSize: '12px', margin: '0 0 8px 0'}
);
 
var readingText = ui.Label(
  'The map shows how each year departs from the 1950-1979 baseline. ' +
  'Blue indicates calling predicted earlier than baseline, red later.',
  {fontSize: '12px', margin: '0 0 12px 0'}
);
 
var instructionText = ui.Label(
  'Click any pixel for its full time series. Use the chart menu to download CSV.',
  {fontSize: '12px', color: '#666', margin: '8px 0 4px 0'}
);
 
var sourcesHeader = ui.Label('Data sources', {
  fontSize: '12px', fontWeight: 'bold', margin: '12px 0 4px 0'
});
 
var sourcesText = ui.Label(
  'Temperature: ERA5-Land Hourly (ECMWF), ~9 km resolution. ' +
  'Range: IUCN Red List, Pseudacris crucifer. ' +
  'Model: Lovett, G.M. (2013) Northeastern Naturalist 20(2):333–340.',
  {fontSize: '12px', margin: '12px 4px 0'}
);
 
// Legend
var legendHeader = ui.Label('Predicted change in calling onset from baseline (days)', {
  fontSize: '12px', fontWeight: 'bold', margin: '0 0 4px 0'
});
 
var legendBar = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
 
var buildLegend = function() {
  var stops = [
    ['#2166ac', '\u2212' + anomalyLimit],
    ['#67a9cf', ''],
    ['#f7f7f7', '0'],
    ['#ef8a62', ''],
    ['#b2182b', + anomalyLimit]
  ];
  stops.forEach(function(stop) {
    legendBar.add(ui.Panel({
      widgets: [ui.Label(stop[1], {
        fontSize: '10px', margin: '2px 0 0 0', textAlign: 'center',
        backgroundColor: 'rgba(0,0,0,0)'
      })],
      style: {backgroundColor: stop[0], width: '60px', padding: '2px'}
    }));
  });
};
 
// Controls
var yearLabel = ui.Label('Year: 1950', {
  fontSize: '14px', fontWeight: 'bold', margin: '12px 0 4px 0'
});
 
var yearSlider = ui.Slider({
  min: 1950, max: 2025, step: 1, value: 1950,
  style: {stretch: 'horizontal', margin: '0 0 12px 0'}
});
 
var chartPanel = ui.Panel({style: {margin: '0'}});
 
// Behaviour
var showYear = function(year) {
  yearLabel.setValue('Year: ' + year);
  var img = ee.Image(anomalyList.get(year - 1950));
  map.layers().set(layerAnomaly, ui.Map.Layer(
    img.clip(geom), anomalyVis, 'Predicted DFC Change ' + year
  ));
};
 
var showSeries = function(coords) {
  chartPanel.clear();
  chartPanel.add(ui.Label('Loading…', {fontSize: '12px', color: '#666'}));
 
  var pt = ee.Geometry.Point([coords.lon, coords.lat]);
  map.layers().set(layerMarker, ui.Map.Layer(pt, {color: 'black'}, 'Selected'));
 
  var chart = ui.Chart.image.series({
    imageCollection: annualDFCCollection.select('predicted_dfc'),
    region: pt,
    reducer: ee.Reducer.first(),
    scale: 9000,
    xProperty: 'system:time_start'
  }).setOptions({
    title: coords.lat.toFixed(2) + '\u00B0N, ' + coords.lon.toFixed(2) + '\u00B0W',
    hAxis: {title: 'Year'},
    vAxis: {title: 'Predicted DFC (day of year)'},
    pointSize: 3,
    lineWidth: 0,
    legend: {position: 'none'},
    trendlines: {0: {showR2: true, visibleInLegend: true, color: '#b2182b'}}
  });
 
  chartPanel.clear();
  chartPanel.add(chart);
};
 
// Trigger changes
yearSlider.onChange(showYear);
map.onClick(showSeries);
 
// Put everything together
buildLegend();
 
mainPanel.add(title);
mainPanel.add(subtitle);
mainPanel.add(methodText);
mainPanel.add(readingText);
mainPanel.add(legendHeader);
mainPanel.add(legendBar);
mainPanel.add(yearLabel);
mainPanel.add(yearSlider);
mainPanel.add(instructionText);
mainPanel.add(chartPanel);
mainPanel.add(sourcesHeader);
mainPanel.add(sourcesText);
 
ui.root.add(mainPanel);
ui.root.add(map);
 
// Set defaults
map.setCenter(-75.5, 43.0, 6);
map.style().set('cursor', 'crosshair');
showYear(1950);
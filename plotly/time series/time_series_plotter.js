const COMMON_LAYOUT = {
  height: 730,
  xaxis: {
    title: 'Year',
  },
  yaxis: {title: 'Billion US$'},
};
const SOURCES = ['SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI',];
const REGIONS = [
  'BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF',
  'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST',
];
const SPECIES_WITH_LEGEND =
    ['0/500', 'CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3',];
const YEARS = [
  1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
  2010, 2011, 2012, 2013, 2014, 2015, 2016,
];
const COLORS = {
  // Species
  CO2: '#ff0066',
  CH4: '#660066',
  BC: '#6600cc',
  SO2: '#0000ff',
  CO: '#006666',
  OC: '#009900',
  N2O: '#ffff00',
  NOx: '#ff3300',
  NH3: '#800026',
  // Regions
  BONA: '#ff0066',
  TENA: '#660066',
  CEAM: '#6600cc',
  NHSA: '#0000ff',
  SHSA: '#006666',
  EURO: '#00004d',
  MIDE: '#009900',
  NHAF: '#ffff00',
  SHAF: '#ff3300',
  BOAS: '#800026',
  TEAS: '#333300',
  SEAS: '663300',
  EQAS: '#999966',
  AUST: '#000000',
  // Sources
  SAVA: '#ff0066',
  BORF: '#660066',
  TEMF: '#0000ff',
  DEFO: '#009900',
  PEAT: '#ff3300',
  AGRI: '#800026',
};
const COLOR_SPECTRUM = [
  [0, '#ffffdd'],
  [.004, '#ffeda0'],
  [.01, '#fed976'],
  [.02, '#feb24c'],
  [.04, '#fd8d3c'],
  [.1, '#fc4e2a'],
  [.2, '#e31a1c'],
  [.4, '#bd0026'],
  [1, '#800026'],
];

const STANDARD_FILL_TRACE = {
  showlegend: false,
  hoverinfo: 'none',
};

const LAST_FILL_TRACE = {
  showlegend: true,
  hoverinfo: 'all',
};

let graphType = 'Area';
let axis = 'Species';
let ENSO = false;

/**
 * Sets the graph type to the selected type.
 * @param {string} graph The selected graph type.
 */
function checkType(graph) {
  graphType = graph;
}

/**
 * Sets the primary axis to the selection.
 * @param {string} selectedAxis The selected primary axis.
 */
function checkAxis(selectedAxis) {
  axis = selectedAxis;
}

/**
 * Sets the ENSO flag, if selected.
 * @param {string} ensoFlag The value of the ENSO flag.
 */
function checkENSO(ensoFlag) {
  ENSO = (ensoFlag == 'ENSO');
}

/**
 * Sums the area underneath each trace for an area/stacked graph. This is so
 * that the area appears to stack cumulatively.
 * @param {!Array<!Object>} traces All the desired traces for a stacked graph.
 * @return {!Array<!Object>} The traces, modified to appear to stack.
 */
function getStackedArea(traces) {
  for (i = 1; i < traces.length; i++) {
    traces[i].y = traces[i].y.map((yearData, yearNumber) =>
        yearData + traces[i - 1].y[yearNumber]);
  }
  return traces;
}

/**
 * Draws an area graph with the given data.
 * @param {!Object<string, string>} primaryList List of primary data sources.
 * @param {string} axis The primary axis in the resultant graph.
 */
function drawStackedGraph(primaryList, axis) {
  // Graph lines, including both filled and unfilled lines.
  const traces = [];

  primaryList.forEach((primaryDatum, primaryIndex) => {
    const primaryName = primaryDatum.name;
    primaryDatum.data.forEach((subDatum, subIndex) => {
      const subName = subDatum.name;
      if (primaryName == 'GLOBAL' || primaryName == 'ALL' || subName == 'GLOBAL') {
        // Don't double our results!
        return;
      }
      let params = STANDARD_FILL_TRACE;
      let name = `${primaryName} ${subName}`;
      if (subIndex == primaryDatum.data.length - 2) {
        params = LAST_FILL_TRACE;
        name = primaryDatum.name;
      }
      const trace = {
        x: YEARS,
        y: subDatum.data,
        marker: {color: COLORS[primaryDatum.name]},
        name: name,
        fill: 'tonexty',
      };
      traces.push(Object.assign(trace, params));
    });
  });

  // Finalize the plot.
  const layout =
      Object.assign(COMMON_LAYOUT, { title: `Annual cost by ${axis}` });
  Plotly.newPlot('myDiv', getStackedArea(traces), layout);
}

/**
 * Draws a line graph with the given data.
 * @param {!Object<string, string>} primaryList List of primary data sources.
 * @param {string} axis The primary axis in the resultant graph.
 */
function drawLineGraph(primaryList, axis) {
  const traces = [];

  primaryList.forEach((primaryData) => {
    if (primaryData.name == 'GLOBAL' || primaryData.name == 'ALL') {
      // Don't double our results!
      return;
    }
    // Vertically sum each sub-matrix into a single object.
    // Result contains 'Global' data for each subject (e.g. "All CO2 emissions")
    const summedDataList = primaryData.data.reduce((summedData, nextData) => {
      if (nextData.name == 'GLOBAL' || nextData.name == 'ALL') {
        return summedData;
      }

      // Add each data point in the next species/source/region to running sum.
      const incrementalSum = summedData.data.map((sum, index) =>
          sum + nextData.data[index]);
      return {
        name: 'ALL',
        data: incrementalSum,
      };
    });

    // Add a line representing all data for this species/source/region.
    traces.push({
      x: YEARS,
      y: summedDataList.data,
      marker: {color: COLORS[primaryData.name]},
      hoverinfo: 'all',
      showlegend: true,
      name: `All ${primaryData.name} emissions`
    });
  });

  // Finalize the plot.
  const layout =
      Object.assign(COMMON_LAYOUT, { title: `Annual cost by ${axis}` });
  Plotly.newPlot('myDiv', traces, layout);
}

/**
 * Draws a heat map with the given data.
 * @param {!Object<string, string>} fileData All the data pulled from the file.
 * @param {string} axis The primary axis in the resultant graph.
 */
function drawHeatMap(fileData, axis) {
  const data = [{
    z: fileData,
    y: SPECIES_WITH_LEGEND,
    type: 'heatmap',
    colorscale: COLOR_SPECTRUM,
    colorbar: {title: 'Billion US$'}
  }];
  const layout = {
    xaxis: {title: axis, side: 'top'},
    yaxis: {title: 'Species'},
    height: 750,
    margin: {t: 150, l: 90, pad: 2},
    pad: 25
  };
  if (axis == 'Regions' || axis == 'Species') {
    data[0].x = REGIONS;
    layout.width = 880;
  } else if (axis == 'Sources') {
    data[0].x = SOURCES;
    layout.width = 530;
  }

  if (ENSO) {
    layout.title = 'SCAR Cost in El Nino 1997-1998';
  } else {
    layout.title = 'Average SCAR Cost 1997-2014';
  }
  Plotly.newPlot('myDiv', data, layout);
}

/**
 * Draws the selected graph type with the given data.
 * @param {!Object<string, string>} fileData All the data pulled from the file.
 * @param {string} axis The primary axis in the resultant graph.
 */
function drawGraph(fileData, axis) {
  const parsedData = JSON.parse(fileData);
  if (graphType == 'Area') {
    drawStackedGraph(parsedData.data, axis);
  } else if (graphType == 'Line') {
    drawLineGraph(parsedData.data, axis);
  } else {
    drawHeatMap(parsedData, axis);
  }
}

/**
 * Loads data from the selected file and draws the designated graph type.
 */
function loadData() {
  const file = fileInput.files[0];
  const textType = /text.*/;

  if (file.type.match(textType)) {
    const reader = new FileReader();
    reader.onload = (event) => drawGraph(reader.result, axis);
    reader.readAsText(file);
  } else {
    alert('Not supported!');
  }
}

/**
 * Set up the window.
 */
window.onload = function() {
  const fileInput = document.getElementById('fileInput');
  fileInput.addEventListener('change', (event) => loadData());
};

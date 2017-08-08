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
  visible: false,
};

const LAST_FILL_TRACE = {
  showlegend: true,
  hoverinfo: 'all',
};

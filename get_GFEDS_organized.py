"""In this example we will calculate annual CO emissions, SCAR, and AQ costs for
the 14 GFED basis regions over 1997-2014. Please adjust the code to calculate
emissions for your own species, region, and time period of interest. Please
first download the GFED4.1s files and the GFED4_Emission_Factors.txt
to your computer and adjust the directory where you placed them below
"""

import csv
import h5py  # Please make sure you have the h5py library
import numpy as np
np.set_printoptions(threshold=np.inf)

months = '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'
"""Strings that are used to represent months in source data."""

regions = [
    'BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF',
    'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST', 'Global'
]
"""BONA: Boreal North America
TENA: Temperate North America
CEAM: Central America
NHSA: Northern Hemisphere South America
SHSA: Southern Hemisphere South America
EURO: Europe
MIDE: Middle East
NHAF: Northern Hemisphere Africa
SHAF: Souther Hemisphere Africa
BOAS: Boreal Asia
TEAS: Temperate Asia
SEAS: Southeast Asia
EQAS: Equatorial Asia
AUST: Australia et al
"""

sources = ['SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI', 'All sources']
"""SAVA: Savanna, grassland, and shrubland fire
BORF: Boreal forest fires
TEMF: Temperate forest fires
DEFO: Tropical deforestation & degradation
PEAT: Peat fires
AGRI: Agricultural waste burning
"""

species_used = ['CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3']
"""The species that we care enough about to plot."""
species_row = 2, 4, 13, 14, 3, 12, 8, 7, 33
"""Rows that the above species are stored on in the source data."""

scar_values = [84, 4589, 270419, 41968, 632, 68299, 36987, 67141, 24906]
"""SCAR scaling factors - eg 1 ton of CO2 = $84 in damages."""

aq_values = [0, 665, 61701, 33009, 239, 50619, 0, 66854, 22432]
"""AQ scaling factors - eg 1 ton of CO2 = $0 in damages."""

regional_scaling = [
    0.09, 3.03, 1.57, 0.37, 0.38, 5.28, 0.28, 1.31, 0.61, 0.45, 7.18, 5.39,
    2.36, 0.15
]
"""Population adjustment factors for AQ scaling."""
global_scaling = 8.53
"""Global population scaling factor."""

data_types = ['emissions', 'SCAR', 'AQ']
units = ['1E12 g', '2007 US$', '2007 US$']
start_year = 1997
end_year = 2014
GRAMS_PER_TON = 1000000
NUM_YEARS = 18
NUM_MONTHS = 12
NUM_SOURCES = 6
NUM_SPECIES = 9
NUM_REGIONS = 15


def read_emissions_factors(directory):
  """Read Emission factors in g species per kg dry matter (DM) burned. Note that
  NOx is as NO. Emission factors based mostly on Akagi et al.
  (2011, http://www.atmos-chem-phys.net/11/4039/2011/acp-11-4039-2011.html),
  but additional data sources were used as well.
  """
  # 41 total species (though we only use 9), 6 sources
  emission_factors = np.zeros((41, 6))

  k = 0
  f = open(directory + '/GFED4_Emission_Factors.txt')
  while 1:
    line = f.readline()
    if line == '':
      break

    if line[0] != '#':
      contents = line.split()
      emission_factors[k, :] = contents[1:]
      k += 1

  f.close()
  return emission_factors


def get_year_file(directory, year):
  """Loads one particular HDF5 file containing dry matter emissions data."""
  string = directory + '/GFED4.1s_' + str(year) + '.hdf5'
  f = h5py.File(string, 'r')
  return f


def load_data(directory, emission_factors):
  """Loads data from HDF5 source dry matter emissions data, emission-species
  proportions, and SCAR/AQ adjustment files into mem.
  """
  # 18 years, 12 months, 6 sources, 9 species, 14 regions
  emissions_data = np.zeros((NUM_YEARS, NUM_MONTHS, NUM_SOURCES, NUM_SPECIES,
                             NUM_REGIONS))

  # 9 species
  for year in range(start_year, end_year + 1):
    print ' '
    print 'Year: ' + str(year)

    f = get_year_file(directory, year)
    basis_regions = f['/ancill/basis_regions'][:]
    grid_area = f['/ancill/grid_cell_area'][:]

    for month in range(NUM_MONTHS):
      print ' '
      print 'Month: ' + str(month)
      # read in Dry Matter emissions
      string = '/emissions/' + months[month] + '/DM'
      dry_matter_emissions = f[string][:]

      for source in range(NUM_SOURCES):
        print 'Source: ' + sources[source]
        # read in the fractional contribution of each source
        string = '/emissions/' + months[month] + '/partitioning/DM_' + sources[source]
        contribution = f[string][:]
        saved_contribution = contribution

        for species_num in range(NUM_SPECIES):
          print 'Species: ' + species_used[species_num]

          # Savannah and Agricultural fires are part of short term carbon cycle.
          if ((sources[source] == 'SAVA' or sources[source] == 'AGRI') and
              species_used[species_num] == 'CO2'):
            contribution = 0
          else:
            contribution = saved_contribution

          source_emissions = np.zeros((720, 1440))
          species_emission_factors = emission_factors[species_row[species_num]]

          source_emissions = dry_matter_emissions * contribution * species_emission_factors[source]

          for region in range(NUM_REGIONS):
            # Don't double count - just multiply by 1.
            if region == NUM_REGIONS - 1:
              mask = np.ones((720, 1440))
            else:
              mask = basis_regions == (region + 1)
            emissions_data[year - start_year, month, source, species_num,
                           region] = np.sum(grid_area * mask * source_emissions)
    print ' '
    print 'YEAR ' + str(year) + ' COMPLETE'

  return emissions_data


def process_scar(value, region_num, species_num):
  """Scales SCAR based on the species and emission region."""
  # AQ effect depends on local population, SCAR should NOT.
  aq_overcompensation = ((
      global_scaling - regional_scaling[region_num]) / global_scaling) * (
          value * aq_values[species_num] / GRAMS_PER_TON)
  return (
      value * scar_values[species_num] / GRAMS_PER_TON) - aq_overcompensation


def process_aq(value, region_num, species_num):
  """Scales air quality costs based on the species and emission region."""
  return (regional_scaling[region_num] / global_scaling) * (
      value * aq_values[species_num] / GRAMS_PER_TON)


def convert_emissions(value):
  """Converts grams to Tg."""
  return value / 1E12


def process_scar_per_ton(value, species_num, carbon_value):
  """Computes the SCAR of a ton of emissions of one species."""
  tons_carbon = carbon_value / GRAMS_PER_TON
  return process_scar(value, species_num) / tons_carbon


def process_aq_per_ton(value, species_num, carbon_value):
  """Computes the Air Quality costs of a ton of emissions of one species."""
  tons_carbon = carbon_value / GRAMS_PER_TON
  return process_aq(value, species_num) / tons_carbon


def plot_all_species_for_year(process_method, emissions_data, year,
                              start_month):
  """Totals species emissions - could be used for a Bokeh bar chart."""
  this_year = year
  # Source (+ all sources) -> species
  all_species_chart = np.zeros((NUM_SPECIES, NUM_SOURCES + 1))
  for species_num in range(NUM_SPECIES):
    totaled_sources = 0
    for source in range(NUM_SOURCES):
      totaled_source = 0
      # Exclude global region
      for region in range(NUM_REGIONS - 1):
        for month in range(NUM_MONTHS):
          this_month = month
          # Deal with when we start partway through a year - such as for ENSO
          if (month + start_month > NUM_MONTHS):
            this_month = ((month + start_month) % NUM_MONTHS) - 1
            this_year = year + 1
          totaled_source += process_method(emissions_data[
              this_year, this_month, source, species_num, region], region,
                                           species_num)
      all_species_chart[species_num, source] = totaled_source
      totaled_sources += totaled_source
    all_species_chart[species_num, NUM_SOURCES] = totaled_sources

  return all_species_chart


def plot_all_regions_for_year(process_method, emissions_data, year,
                              start_month):
  """Totals regional emissions - could be used for a Bokeh bar chart."""
  this_year = year
  # source (+ all sources) -> regions (except global)
  all_regions_chart = np.zeros((NUM_REGIONS - 1, NUM_SOURCES + 1))
  # exclude global region
  for region in range(NUM_REGIONS - 1):
    totaled_sources = 0
    for source in range(NUM_SOURCES):
      totaled_source = 0
      for species_num in range(NUM_SPECIES):
        for month in range(NUM_MONTHS):
          this_month = month
          # dealing with when we start partway through a year - such as for ENSO
          if (month + start_month > NUM_MONTHS):
            this_month = ((month + start_month) % NUM_MONTHS) - 1
            this_year = year + 1
          totaled_source += process_method(emissions_data[
              this_year, this_month, source, species_num, region], region,
                                           species_num)
      all_regions_chart[region, source] = totaled_source
      totaled_sources += totaled_source
    all_regions_chart[region, NUM_SOURCES] = totaled_sources

  return all_regions_chart


def plot_regions_species_for_year(process_method, emissions_data, year,
                                  start_month):
  """Totals emissions for one region, in one year."""
  this_year = year
  # source (+ all sources) -> regions (except global)
  combined_chart = np.zeros((NUM_SPECIES, NUM_REGIONS))
  # exclude global region
  for species_num in range(NUM_SPECIES):
    totaled_species = 0
    for region in range(NUM_REGIONS - 1):
      totaled_region = 0
      for source in range(NUM_SOURCES):
        for month in range(NUM_MONTHS):
          this_month = month
          # dealing with when we start partway through a year - such as for ENSO
          if month + start_month > NUM_MONTHS:
            this_month = ((month + start_month) % NUM_MONTHS) - 1
            this_year = year + 1
          totaled_region += process_method(emissions_data[
              this_year, this_month, source, species_num, region], region,
                                           species_num)
      combined_chart[species_num, region] = totaled_region
      totaled_species += totaled_region
    combined_chart[species_num, NUM_REGIONS - 1] = totaled_species

  return combined_chart


def plot_time_series_for_sources(emissions_data, process_method):
  """Computes emissions for one source over entire time series."""
  for source in range(NUM_SOURCES):
    if sources[source] != 'SAVA':
      all_years_chart = np.zeros((NUM_YEARS, NUM_REGIONS - 1))
      for year in range(NUM_YEARS):
        for region in range(NUM_REGIONS - 1):
          totaled_region = 0
          totaled_carbon = 0
          for species_num in range(NUM_SPECIES):
            for month in range(NUM_MONTHS):
              totaled_carbon += emissions_data[year, month, source, 0,
                                               region] / GRAMS_PER_TON
              totaled_region += process_method(
                  emissions_data[year, month, source, species_num, region],
                  region, species_num)
          if totaled_carbon > 0:
            all_years_chart[year, region] = totaled_region / totaled_carbon
          return all_years_chart


def write_species(emissions_data, data_type, process_method, species_num,
                  units):
  """Write data for a single species to CSV file."""
  writer = csv.writer(
      open('./plots/tables/' + data_type + '/' + species_used[species_num] +
           '.csv', 'w'))
  writer.writerow(['units - ' + units])
  writer.writerow([species_used[species_num]])
  for year in range(NUM_YEARS):
    writer.writerow([])
    writer.writerow([start_year + year])
    region_header = [
        '', 'BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF',
        'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST', 'Global'
    ]
    writer.writerow(region_header)
    totaled_regions = ['total', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for source in range(NUM_SOURCES):
      region_list = [
          sources[source], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      ]
      totaled_source = 0
      for region in range(NUM_REGIONS - 1):
        totaled_region = 0
        for month in range(NUM_MONTHS):
          totaled_region += process_method(
              emissions_data[year, month, source, species_num, region], region,
              species_num)
        region_list[region + 1] = totaled_region
        totaled_regions[region + 1] += totaled_region
        totaled_source += totaled_region
      region_list[NUM_REGIONS] = totaled_source
      totaled_regions[NUM_REGIONS] += totaled_source
      writer.writerow(region_list)
    writer.writerow(totaled_regions)


def write_data(emissions_data):
  """Write all data to CSV files."""
  for species in range(NUM_SPECIES):
    write_species(emissions_data, 'air_quality', process_aq, species,
                  '2007 US$')
    write_species(emissions_data, 'emissions', convert_emissions, species,
                  '1E12 g')
    write_species(emissions_data, 'SCAR', process_scar, species, '2007 US$')


"""Print time series and heatmap data to the console"""


def plot_all_years(identifier, emissions_data):
  emissions_time_series = []
  scar_time_series = []
  air_quality_time_series = []
  scar_regions_species = []
  aq_regions_species = []
  # each calendar year -- append to a time series array
  for year in range(NUM_YEARS):
    emissions_time_series.append(
        plot_all_species_for_year('emissions', convert_emissions,
                                  emissions_data, year, 0,
                                  str(year + start_year)))
    scar_time_series.append(
        plot_all_species_for_year('SCAR', process_scar, emissions_data, year, 0,
                                  str(year + start_year)))
    air_quality_time_series.append(
        plot_all_species_for_year('air_quality', process_aq, emissions_data,
                                  year, 0, str(year + start_year)))
    scar_regions_species.append(
        plot_regions_species_for_year('SCAR', process_scar, emissions_data,
                                      year, 0, str(year + start_year)))
    aq_regions_species.append(
        plot_regions_species_for_year('air_quality', process_aq, emissions_data,
                                      year, 0, str(year + start_year)))
  # heatmaps for 97-98
  ENSO_SCAR_source_species = plot_all_species_for_year(
      'SCAR', process_scar, emissions_data, 1997 - start_year, 7, 'El Nino')
  ENSO_AQ_source_species = plot_all_species_for_year(
      'air_quality', process_aq, emissions_data, 1997 - start_year, 7,
      'El Nino')
  ENSO_SCAR_regions_species = plot_regions_species_for_year(
      'SCAR', process_scar, emissions_data, 1997 - start_year, 7, 'El Nino')
  ENSO_AQ_regions_species = plot_regions_species_for_year(
      'air_quality', process_aq, emissions_data, 1997 - start_year, 7,
      'El Nino')
  print 'ENSO SCAR SOURCE SPECIES'
  print ENSO_SCAR_source_species
  print 'ENSO AQ SOURCE SPECIES'
  print ENSO_AQ_source_species
  print 'ENSO SCAR REGIONS SPECIES'
  print ENSO_SCAR_regions_species
  print 'ENSO AQ REGIONS SPECIES'
  print ENSO_AQ_regions_species
  # heatmaps for the average of all years
  print 'avg'
  print '\nAVERAGES\n'
  print 'EMISSIONS AVERAGES\n'
  print np.divide(np.sum(emissions_time_series, axis=0), NUM_YEARS)
  print '\nSCAR AVERAGES\n'
  print np.divide(np.sum(scar_time_series, axis=0), NUM_YEARS)
  print '\nAIR QUALITY\n'
  print np.divide(np.sum(air_quality_time_series, axis=0), NUM_YEARS)
  print '\nSCAR AVERAGES -- SPECIES AND REGIONS\n'
  print np.divide(np.sum(scar_regions_species, axis=0), NUM_YEARS)
  print '\nAIR QUALITY AVERAGES -- SPECIES AND REGIONS\n'
  print np.divide(np.sum(aq_regions_species, axis=0), NUM_YEARS)
  # time series that show all years in one chart
  print '\nSCAR SOURCES\n'
  print np.transpose(scar_time_series, (2, 1, 0))
  print '\nSCAR REGIONS\n'
  print np.transpose(scar_regions_species, (2, 1, 0))
  print '\nSCAR SPECIES\n'
  print np.transpose(scar_regions_species, (1, 2, 0))
  print '\nAIR QUALITY SOURCES\n'
  print np.transpose(air_quality_time_series, (2, 1, 0))
  print '\nAIR QUALITY REGIONS\n'
  print np.transpose(aq_regions_species, (2, 1, 0))
  print '\nAIR QUALITY SPECIES\n'
  print np.transpose(aq_regions_species, (1, 2, 0))


def plot_data(emissions_data):
  """Plot all data using Bokeh."""
  plot_all_years('species', emissions_data)


def calculate_emissions():
  """Calculate, plot, and write emissions and cost data."""
  directory = '.'

  emission_factors = read_emissions_factors(directory)

  print 'loading...'
  emissions_data = load_data(directory, emission_factors)

  print ''
  print 'plotting...'
  plot_data(emissions_data)
  write_data(emissions_data)


# compare this to http://www.falw.vu/~gwerf/GFED/GFED4/tables/GFED4.1s_CO.txt
if __name__ == '__main__':
  calculate_emissions()

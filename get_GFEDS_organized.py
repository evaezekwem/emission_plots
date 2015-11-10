from plot import Plotter
import csv
import numpy as np
import h5py # if this creates an error please make sure you have the h5py library

months       = '01','02','03','04','05','06','07','08','09','10','11','12'
regions      = ['BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST', 'Global'];
sources      = ['SAVA','BORF','TEMF','DEFO','PEAT','AGRI', 'All sources'];
species_used = ['CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3']
species_row  = 2, 4, 13, 14, 3, 12, 8, 7, 33
scar_values = [84, 4589, 270419, 41968, 632, 68299, 36987, 67141, 24906];
aq_values = [0, 665, 61701, 33009, 253, 50619, 0, 66854, 22432];
data_types = ["emissions", "SCAR", "AQ"];
units = ["1E12 g", "2007 US$", "2007 US$"];
start_year = 1997
end_year   = 2014
GRAMS_PER_TON = 907185;
NUM_YEARS = 18;
NUM_MONTHS = 12;
NUM_SOURCES = 6;
NUM_SPECIES = 9;
NUM_REGIONS = 15;

def read_emissions_factors(directory):
    """
    Read in emission factors
    """
    EFs     = np.zeros((41, 6)) # 41 species, 6 sources
    
    k = 0
    f = open(directory+'/GFED4_Emission_Factors.txt')
    while 1:
        line = f.readline()
        if line == "":
            break
            
        if line[0] != '#':
            contents = line.split()
            EFs[k,:] = contents[1:]
            k += 1
                    
    f.close()
    return EFs;

def load_data(directory, EFs):
    # 18 years, 12 months, 6 sources, 9 species, 14 regions
    emissions_data = np.zeros((NUM_YEARS, NUM_MONTHS, NUM_SOURCES, NUM_SPECIES, NUM_REGIONS));
    
    # 9 species
    for year in range(start_year, end_year):
        print "year: " + str(year);

        f = get_year_file(directory, year);
        basis_regions = f['/ancill/basis_regions'][:]
        grid_area     = f['/ancill/grid_cell_area'][:]
        
        for month in range(NUM_MONTHS):
            print "month: " + str(month);
            # read in DM emissions
            string = '/emissions/'+months[month]+'/DM'
            DM_emissions = f[string][:]
            
            for source in range(NUM_SOURCES):
                # read in the fractional contribution of each source
                string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
                contribution = f[string][:]
            
                for species_num in range(NUM_SPECIES):
                    print " "
                    print "Species: " + species_used[species_num]
                    
                    if(sources[source] == 'SAVA' and species_used[species_num] == 'CO2'):
                        contribution = 0;
                    
                    source_emissions = np.zeros((720, 1440))
                    EF_species = EFs[species_row[species_num]];
                    
                    source_emissions = DM_emissions * contribution * EF_species[source];
                    
                    for region in range(NUM_REGIONS):
                        #global
                        if region == NUM_REGIONS - 1:
                            mask = np.ones((720, 1440))
                        else:
                            mask = basis_regions == (region + 1) 
                        emissions_data[year - start_year, month, source, species_num, region] = np.sum(grid_area * mask * emissions_data);
                        
    return emissions_data;
    
def finalize_and_plot(table, plotter_method, year_label):
    # convert to $ value
    scar_table = table * scar_values[species_num] / GRAMS_PER_TON;
    aq_table = table * aq_values[species_num] / GRAMS_PER_TON;
    # convert to Tg CO 
    emissions_table = table / 1E12;
    plot_species_total(plotter, year_label + "_all_species", "SCAR", scar_table);
    plot_species_total(plotter, year_label + "_all_species", "air_quality", aq_table);
    plot_species_total(plotter, year_label + "_all_species", "emissions", emissions_table);
    
def plot_all_species_for_year(start_year, start_month, year_label):
    this_year = start_year;
    #species -> source (+ all sources)
    all_species_chart = np.zeros(NUM_SPECIES, NUM_SOURCES + 1);
    for species_num in range(NUM_SPECIES):
        totaled_sources = 0;
        for source in range (NUM_SOURCES):
            totaled_source = 0;
            for region in range(NUM_REGIONS):
                for month in range(NUM_MONTHS):
                    this_month = month
                    # dealing with when we start partway through a year -- such as for ENSO
                    if(month + start_month > NUM_MONTHS):
                        this_month = month % NUM_MONTHS;
                        this_year = start_year + 1;
                    totaled_source += emissions_data[this_year, this_month, source, species_num, region];
            all_species_chart[species_num, source] = totaled_source;
            totaled_sources += totaled_source;
        all_species_chart[species_num, NUM_SOURCES] = totaled_sources;
        
    finalize_and_plot(all_species_chart, plotter.plot_species_total, year_label);
    
def plot_all_regions_for_year(start_year, start_month, year_label):
    this_year = start_year;
    #species -> source (+ all sources)
    all_regions_chart = np.zeros(NUM_REGIONS, NUM_SOURCES + 1);
    for region in range(NUM_REGIONS):
        totaled_sources = 0;
        for source in range (NUM_SOURCES):
            totaled_source = 0;
            for species_num in range(NUM_SPECIES):
                for month in range(NUM_MONTHS):
                    this_month = month
                    # dealing with when we start partway through a year -- such as for ENSO
                    if(month + start_month > NUM_MONTHS):
                        this_month = month % NUM_MONTHS;
                        this_year = start_year + 1;
                    totaled_source += emissions_data[this_year, this_month, source, species_num, region];
            all_regions_chart[species_num, source] = totaled_source;
            totaled_sources += totaled_source;
        all_regions_chart[species_num, NUM_SOURCES] = totaled_sources;
    finalize_and_plot(all_regions_chart, plotter.plot_regions_total, year_label);
    
def plot_all_years(emissions_data, plot_method):
    # each calendar year
    for year in range(NUM_YEARS):
        plot_method(year, 0, year+start_year);
    # ENSO years -- july to june
    plot_method(1997, 7, "97-98_El_Nino");
    plot_method(1997, 7, "98-99_La_Nina");
            
def plot_data(emissions_data):
    plotter = Plotter();
    
    # each species for each year
    #for year in range(18):
    #    for species_num in range(9):
            
    #        plotter.plot_species() 
    
    # total for all species for each year
    plot_all_years(emissions_data, plot_all_species_for_year);
    
    # total for all regions for each year
    plot_all_years(emissions_data, plot_all_regions_for_year);
    
    # time series
        

def calculate_emissions():
    
    # in this example we will calculate annual CO emissions for the 14 GFED 
    # basisregions over 1997-2014. Please adjust the code to calculate emissions
    # for your own specie, region, and time period of interest. Please
    # first download the GFED4.1s files and the GFED4_Emission_Factors.txt
    # to your computer and adjust the directory where you placed them below
    directory    = '.'

    EFs = read_emissions_factors(directory);
    
    print "loading...";
    emissions_data = load_data(directory, EFs);
    
    print "";
    print "plotting...";
    plot_data(emissions_data);

    


# please compare this to http://www.falw.vu/~gwerf/GFED/GFED4/tables/GFED4.1s_CO.txt
if __name__ == '__main__':
        calculate_emissions()

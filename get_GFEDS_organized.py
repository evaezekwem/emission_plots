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
regional_scaling = [0.09, 3.03, 1.57, 0.37, 0.38, 5.28, 0.28, 1.31, 0.61, 0.45, 7.18, 5.39, 2.36, 0.15];
global_scaling = 8.53;
data_types = ["emissions", "SCAR", "AQ"];
units = ["1E12 g", "2007 US$", "2007 US$"];
start_year = 1997
end_year   = 2014
GRAMS_PER_TON = 1000000;
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

def get_year_file(directory, year):
    string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
    f = h5py.File(string, 'r')
    return f

def load_data(directory, EFs):
    # 18 years, 12 months, 6 sources, 9 species, 14 regions
    emissions_data = np.zeros((NUM_YEARS, NUM_MONTHS, NUM_SOURCES, NUM_SPECIES, NUM_REGIONS));
    
    # 9 species
    for year in range(start_year, end_year + 1):
        print " ";
        print "Year: " + str(year);

        f = get_year_file(directory, year);
        basis_regions = f['/ancill/basis_regions'][:]
        grid_area     = f['/ancill/grid_cell_area'][:]
        
        for month in range(NUM_MONTHS):
            print " ";
            print "Month: " + str(month);
            # read in DM emissions
            string = '/emissions/'+months[month]+'/DM'
            DM_emissions = f[string][:]
            
            for source in range(NUM_SOURCES):
                print "Source: " + sources[source];
                # read in the fractional contribution of each source
                string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
                contribution = f[string][:]
                saved_contribution = contribution;

                for species_num in range(NUM_SPECIES):
                    print "Species: " + species_used[species_num]
                    
                    if((sources[source] == 'SAVA' or sources[source] == 'AGRI') and species_used[species_num] == 'CO2'):
                        contribution = 0;
                    else:
                        contribution = saved_contribution;
                    
                    source_emissions = np.zeros((720, 1440))
                    EF_species = EFs[species_row[species_num]];
                    
                    source_emissions = DM_emissions * contribution * EF_species[source];
                    
                    for region in range(NUM_REGIONS):
                        #global
                        if region == NUM_REGIONS - 1:
                            mask = np.ones((720, 1440))
                        else:
                            mask = basis_regions == (region + 1) 
                        emissions_data[year - start_year, month, source, species_num, region] = np.sum(grid_area * mask * source_emissions);
        print " ";
        print "YEAR " + str(year) + " COMPLETE";
                        
    return emissions_data;

# convert to $ value  
def process_scar(value, region_num, species_num):
    aq_overcompensation = ((global_scaling-regional_scaling[region_num])/global_scaling) * (value *aq_values[species_num] / GRAMS_PER_TON);
    return (value * scar_values[species_num] / GRAMS_PER_TON) - aq_overcompensation;
    
def process_aq(value, region_num, species_num):
    return (regional_scaling[region_num]/global_scaling) * (value * aq_values[species_num] / GRAMS_PER_TON);

# convert to Tg CO  
def process_emissions(value, region_num, species_num):
    return value / 1E12;
    
def process_scar_per_ton(value, species_num, carbon_value):
    tons_carbon = carbon_value / GRAMS_PER_TON;
    return process_scar(value, species_num) / tons_carbon;
    
def process_aq_per_ton(value, species_num, carbon_value):
    tons_carbon = carbon_value / GRAMS_PER_TON;
    return process_aq(value, species_num) / tons_carbon;
    
def plot_all_species_for_year(plotter, metric, process_method, emissions_data, year, start_month, year_label):
    this_year = year;
    #source (+ all sources) -> species
    all_species_chart = np.zeros((NUM_SOURCES + 1, NUM_SPECIES));
    for species_num in range(NUM_SPECIES):
        totaled_sources = 0;
        for source in range (NUM_SOURCES):
            totaled_source = 0;
            #exclude global region
            for region in range(NUM_REGIONS - 1):
                for month in range(NUM_MONTHS):
                    this_month = month
                    # dealing with when we start partway through a year -- such as for ENSO
                    if(month + start_month > NUM_MONTHS):
                        this_month = ((month + start_month) % NUM_MONTHS) - 1;
                        this_year = year + 1;
                    totaled_source += process_method(emissions_data[this_year, this_month, source, species_num, region], region, species_num);
            all_species_chart[source, species_num] = totaled_source;
            totaled_sources += totaled_source;
        all_species_chart[NUM_SOURCES, species_num] = totaled_sources;
        
    plotter.plot_species_total(year_label + "_all_species", metric, all_species_chart);
    return all_species_chart;
    
def plot_all_regions_for_year(plotter, metric, process_method, emissions_data, year, start_month, year_label):
    this_year = year;
    #source (+ all sources) -> regions (except global)
    all_regions_chart = np.zeros((NUM_REGIONS - 1, NUM_SOURCES + 1));
    #exclude global region
    for region in range(NUM_REGIONS - 1):
        totaled_sources = 0;
        for source in range (NUM_SOURCES):
            totaled_source = 0;
            for species_num in range(NUM_SPECIES):
                for month in range(NUM_MONTHS):
                    this_month = month;
                    # dealing with when we start partway through a year -- such as for ENSO
                    if(month + start_month > NUM_MONTHS):
                        this_month = ((month + start_month) % NUM_MONTHS) - 1;
                        this_year = year + 1;
                    totaled_source += process_method(emissions_data[this_year, this_month, source, species_num, region], region, species_num);
            all_regions_chart[region, source] = totaled_source;
            totaled_sources += totaled_source;
        all_regions_chart[region, NUM_SOURCES] = totaled_sources;
        
    plotter.plot_regions_total(year_label + "_all_regions", metric, all_regions_chart);
    return all_regions_chart;
    
def plot_regions_species_for_year(plotter, metric, process_method, emissions_data, year, start_month, year_label):
    this_year = year;
    #source (+ all sources) -> regions (except global)
    combined_chart = np.zeros((NUM_SPECIES, NUM_REGIONS - 1));
    #exclude global region
    for region in range(NUM_REGIONS - 1):
        for species_num in range(NUM_SPECIES):
            totaled_species = 0;
            for source in range (NUM_SOURCES):
                for month in range(NUM_MONTHS):
                    this_month = month;
                    # dealing with when we start partway through a year -- such as for ENSO
                    if(month + start_month > NUM_MONTHS):
                        this_month = ((month + start_month) % NUM_MONTHS) - 1;
                        this_year = year + 1;
                    totaled_species += process_method(emissions_data[this_year, this_month, source, species_num, region], region, species_num);
            combined_chart[species_num, region] = totaled_species;
        
    #plotter.plot_regions_total(year_label + "_species_regions", metric, all_regions_chart);
    return combined_chart;
    
def plot_all_years(plotter, identifier, emissions_data, plot_method, time_series_method):
    emissions_time_series = [];
    scar_time_series = [];
    air_quality_time_series = [];
    scar_regions_species = [];
    aq_regions_species = [];
    # each calendar year -- append to a time series array
    for year in range(NUM_YEARS):
        emissions_time_series.append(plot_method(plotter, "emissions", process_emissions, emissions_data, year, 0, str(year+start_year)));
        scar_time_series.append(plot_method(plotter, "SCAR", process_scar, emissions_data, year, 0, str(year+start_year)));
        air_quality_time_series.append(plot_method(plotter, "air_quality", process_aq, emissions_data, year, 0, str(year+start_year)));
        scar_regions_species.append(plot_regions_species_for_year(plotter, "SCAR", process_scar, emissions_data, year, 0, str(year+start_year)));
        aq_regions_species.append(plot_regions_species_for_year(plotter, "air_quality", process_aq, emissions_data, year, 0, str(year+start_year)))
    # heatmaps for 97-98
    ENSO_SCAR_source_species = plot_method(plotter, "SCAR", process_scar, emissions_data, 1997 - start_year, 7, "El Nino");
    ENSO_AQ_source_species = plot_method(plotter, "air_quality", process_aq, emissions_data, 1997 - start_year, 7, "El Nino");
    ENSO_SCAR_regions_species = plot_regions_species_for_year(plotter, "SCAR", process_scar, emissions_data, 1997 - start_year, 7, "El Nino");
    ENSO_AQ_regions_species = plot_regions_species_for_year(plotter, "air_quality", process_aq, emissions_data, 1997 - start_year, 7, "El Nino");
    print "ENSO SCAR SOURCE SPECIES";
    print ENSO_SCAR_source_species;
    print "ENSO AQ SOURCE SPECIES";
    print ENSO_AQ_source_species;
    print "ENSO SCAR REGIONS SPECIES";
    print ENSO_SCAR_regions_species;
    print "ENSO AQ REGIONS SPECIES";
    print ENSO_AQ_regions_species;
    plotter.plot_heatmap(identifier, False, "SCAR", ENSO_SCAR_source_species);
    plotter.plot_heatmap(identifier, False, "air_quality", ENSO_AQ_source_species);
    plotter.plot_heatmap(identifier, False, "SCAR", ENSO_SCAR_regions_species);
    plotter.plot_heatmap(identifier, False, "air_quality", ENSO_AQ_regions_species);
    # heatmaps for the average of all years
    print "avg"
    print ("\nAVERAGES\n");
    print ("EMISSIONS AVERAGES\n");
    print (np.divide(np.sum(emissions_time_series, axis=0), NUM_YEARS));
    print ("\nSCAR AVERAGES\n");
    print (np.divide(np.sum(scar_time_series, axis=0), NUM_YEARS));
    print ("\nAIR QUALITY\n");
    print (np.divide(np.sum(air_quality_time_series, axis=0), NUM_YEARS));
    print ("\nSCAR AVERAGES -- SPECIES AND REGIONS\n");
    print (np.divide(np.sum(scar_regions_species, axis=0), NUM_YEARS));
    print ("\nAIR QUALITY AVERAGES -- SPECIES AND REGIONS\n");
    print (np.divide(np.sum(aq_regions_species, axis=0), NUM_YEARS));
    plotter.plot_heatmap(identifier, False, "emissions", np.divide(np.sum(emissions_time_series, axis=0), NUM_YEARS));
    plotter.plot_heatmap(identifier, False, "SCAR", np.divide(np.sum(scar_time_series, axis=0), NUM_YEARS));
    plotter.plot_heatmap(identifier, False, "air_quality", np.divide(np.sum(air_quality_time_series, axis=0), NUM_YEARS));
    plotter.plot_heatmap("both", False, "SCAR", np.divide(np.sum(scar_regions_species, axis=0), NUM_YEARS));
    plotter.plot_heatmap("both", False, "air_quality", np.divide(np.sum(aq_regions_species, axis=0), NUM_YEARS));
    # time series that show all years in one chart
    time_series_method(identifier + "_emissions_time_series", "emissions", np.concatenate(emissions_time_series, axis=0));
    time_series_method(identifier + "_SCAR_time_series", "SCAR", np.concatenate(scar_time_series, axis=0));
    time_series_method(identifier + "_air_quality_time_series", "air_quality", np.concatenate(air_quality_time_series, axis=0));
    
    #ENSO years -- july to june
    plot_method(plotter, "emissions", process_emissions, emissions_data, 1997 - start_year, 7, "97-98_El_Nino");
    plot_method(plotter, "SCAR", process_scar, emissions_data, 1997 - start_year, 7, "97-98_El_Nino");
    plot_method(plotter, "air_quality", process_aq, emissions_data, 1997 - start_year, 7, "97-98_El_Nino");
    plot_method(plotter, "emissions", process_emissions, emissions_data, 1998 - start_year, 7, "98-99_La_Nina");
    plot_method(plotter, "SCAR", process_scar, emissions_data, 1998 - start_year, 7, "98-99_La_Nina");
    plot_method(plotter, "air_quality", process_aq, emissions_data, 1998 - start_year, 7, "98-99_La_Nina");
    
# per ton carbon
def plot_time_series_for_sources(plotter, emissions_data, process_method, metric):
    all_sources_all_years_chart = np.zeros((NUM_YEARS, NUM_REGIONS - 1));
    for source in range (NUM_SOURCES):
        if(sources[source] != "SAVA"):
            all_years_chart = np.zeros((NUM_YEARS, NUM_REGIONS - 1));
            for year in range(NUM_YEARS):
                for region in range(NUM_REGIONS - 1):
                    totaled_region = 0;
                    totaled_carbon = 0;
                    #tons_of_carbon = 0;
                    for species_num in range(NUM_SPECIES):
                        for month in range(NUM_MONTHS):
                            totaled_carbon += emissions_data[year, month, source, 0, region] / GRAMS_PER_TON;
                            totaled_region += process_method(emissions_data[year, month, source, species_num, region], region, species_num);
                    if(totaled_carbon > 0):
                        all_years_chart[year, region] = totaled_region / totaled_carbon;
                    #all_sources_all_years_chart[year,region] += totaled_region;
            plotter.plot_regions_source_time_series(sources[source] + "_per_ton_carbon", metric, all_years_chart);
    plotter.plot_regions_source_time_series("All_sources_per_ton_carbon", metric, all_sources_all_years_chart);
            #return all_years_chart;

def write_species(emissions_data, data_type, process_method, species_num, units):
    writer = csv.writer(open("./plots/tables/" + data_type + "/" + species_used[species_num] + ".csv", "w"));
    writer.writerow(["units - " + units]);
    writer.writerow([species_used[species_num]]);
    for year in range(NUM_YEARS):
        writer.writerow([]);
        writer.writerow([start_year+year]);
        region_header = ['', 'BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST', 'Global'];
        writer.writerow(region_header);
        totaled_regions = ["total", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for source in range(NUM_SOURCES):
            region_list = [sources[source], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            totaled_source = 0;
            for region in range(NUM_REGIONS - 1):
                totaled_region = 0;
                for month in range(NUM_MONTHS):
                    totaled_region += process_method(emissions_data[year, month, source, species_num, region], region, species_num);
                region_list[region+1] = totaled_region;
                totaled_regions[region+1] += totaled_region;
                totaled_source += totaled_region;
            region_list[NUM_REGIONS] = totaled_source;
            totaled_regions[NUM_REGIONS] += totaled_source;
            writer.writerow(region_list);
        writer.writerow(totaled_regions);
        
def write_data(emissions_data):
    for species in range(NUM_SPECIES):
        write_species(emissions_data, "air_quality", process_aq, species, "2007 US$");
        write_species(emissions_data, "emissions", process_emissions, species, "1E12 g");
        write_species(emissions_data, "SCAR", process_scar, species, "2007 US$");
            
def plot_data(emissions_data):
    plotter = Plotter();
    
    # each species for each year
    #for year in range(18):
    #    for species_num in range(9):
            
    
    # total for all species for each year
    plot_all_years(plotter, "species", emissions_data, plot_all_species_for_year, plotter.plot_species_time_series);
    
    # total for all regions for each year
    #plot_all_years(plotter, "regions", emissions_data, plot_all_regions_for_year, plotter.plot_regions_time_series);
    
    #        plotter.plot_species() 
    # time series for each source
    # plot_time_series_for_sources(plotter, emissions_data, process_emissions, "emissions");
    plot_time_series_for_sources(plotter, emissions_data, process_scar, "SCAR");
    plot_time_series_for_sources(plotter, emissions_data, process_aq, "air_quality");
   
        

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
    write_data(emissions_data);

# please compare this to http://www.falw.vu/~gwerf/GFED/GFED4/tables/GFED4.1s_CO.txt
if __name__ == '__main__':
        calculate_emissions()

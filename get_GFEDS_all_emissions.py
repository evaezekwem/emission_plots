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

def setup_writer(data_type, species, units):
    writer = csv.writer(open("./plots/tables/" + data_type + "/" + species + '_' + data_type + '.csv', 'w'));
    writer.writerow(["units - " + units]);
    return writer;

#mutator method to load sums of emissions data over a region into the result table
def sum_regions(result_table, grid_area, emissions_data, basis_regions, source):
    for region in range(15):
        if region == 14:
            mask = np.ones((720, 1440))
        else:
            mask = basis_regions == (region + 1)            
        result_table[source, region] = np.sum(grid_area * mask * emissions_data);     

#plots data for a species' impact and writes it to a csv file
def plot_and_write(plotter, writer, result_table, data_type, identifier):
    plotter.plot_regions(identifier, data_type, result_table);
    for source in range(7):
        source_list = result_table[source, :].tolist();
        source_list.insert(0, sources[source]);
        writer.writerow(source_list);
    writer.writerow([]);
    
def get_year_file(directory, year):
    string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
    f = h5py.File(string, 'r')
    return f

#calculates information for a species' impact in a year   
def calculate_species_for_year(directory, species_num, EF_species, year, start_month):
    emissions_table = np.zeros((7, 15)) # source, region
    source_emissions = np.zeros(7);
    f = get_year_file(directory, year);
    basis_regions = f['/ancill/basis_regions'][:]
    grid_area     = f['/ancill/grid_cell_area'][:]
    
    total_emissions = np.zeros((720, 1440));
    for source in range(6):
        source_emissions = np.zeros((720, 1440))
        for month in range(12):
            #since we're not always doing january-december, may have to load a new year
            if(month == 12):
                f = get_year_file(directory, year+1);
            new_month = (month+start_month) % 12;
            # read in DM emissions
            string = '/emissions/'+months[new_month]+'/DM'
            DM_emissions = f[string][:]
            # read in the fractional contribution of each source
            string = '/emissions/'+months[new_month]+'/partitioning/DM_'+sources[source]
            contribution = f[string][:]
            if(sources[source] == 'SAVA' and species_used[species_num] == 'CO2'):
                contribution = 0;
            source_emissions += DM_emissions * contribution * EF_species[source];
        # calculate CO emissions as the product of DM emissions (kg DM per 
        # m2 per month), the fraction the specific source contributes to 
        # this (unitless), and the emission factor (g CO per kg DM burned)
        total_emissions += source_emissions;
        sum_regions(emissions_table, grid_area, source_emissions, basis_regions, source);
    
    # fill table with total values for the globe (row 15) or basisregion (1-14)
    sum_regions(emissions_table, grid_area, total_emissions, basis_regions, 6);
    print "\t" + str(year) + " done";
    return emissions_table;   

#post-processes the data calculated for a species impact in a year and writes to all the different data files      
def plot_and_write_table(emissions_table, writers, plotter, species_num, EF_species, identifier):
    for writer in writers:
        writer.writerow([identifier]);
        writer.writerow([""] + regions);

    # convert to $ value
    scar_table = emissions_table * scar_values[species_num] / GRAMS_PER_TON;
    aq_table = emissions_table * aq_values[species_num] / GRAMS_PER_TON;
    # convert to Tg CO 
    final_emissions_table = emissions_table / 1E12;
    tables = [final_emissions_table, scar_table, aq_table];
    for data_type in range(3):
        plot_and_write(plotter, writers[data_type], tables[data_type], data_types[data_type], identifier);
    print emissions_table

def plot_regions_table(emissions_table, plotter, species_num, EF_species, identifier):
        # convert to $ value
    scar_table = emissions_table * scar_values[species_num] / GRAMS_PER_TON;
    aq_table = emissions_table * aq_values[species_num] / GRAMS_PER_TON;
    # convert to Tg CO 
    final_emissions_table = emissions_table / 1E12;
    tables = [final_emissions_table, scar_table, aq_table];
    for data_type in range(3):
        plotter.plot_regions(identifier, data_types[data_type], tables[data_type]);
    print emissions_table

def plot_species_table(emissions_table, plotter, species_num, EF_species, identifier):
        # convert to $ value
    scar_table = emissions_table * scar_values[species_num] / GRAMS_PER_TON;
    aq_table = emissions_table * aq_values[species_num] / GRAMS_PER_TON;
    # convert to Tg CO 
    final_emissions_table = emissions_table / 1E12;
    tables = [final_emissions_table, scar_table, aq_table];
    for data_type in range(3):
        plotter.plot_species(identifier, data_types[data_type], tables[data_type]);
    print emissions_table

def calculate_emissions():
    
    # in this example we will calculate annual CO emissions for the 14 GFED 
    # basisregions over 1997-2014. Please adjust the code to calculate emissions
    # for your own specie, region, and time period of interest. Please
    # first download the GFED4.1s files and the GFED4_Emission_Factors.txt
    # to your computer and adjust the directory where you placed them below
    directory    = '.'


    """
    Read in emission factors
    """
    species = [] # names of the different gas and aerosol species
    EFs     = np.zeros((41, 6)) # 41 species, 6 sources

    k = 0
    f = open(directory+'/GFED4_Emission_Factors.txt')
    while 1:
        line = f.readline()
        if line == "":
            break
            
        if line[0] != '#':
            contents = line.split()
            species.append(contents[0])
            EFs[k,:] = contents[1:]
            k += 1
                    
    f.close()

    plotter = Plotter();
    regional_totals = [np.zeros((7, 15))] * 18;
    species_totals = [np.zeros((7, 9))] * 18;
    for species_num in range(9):
        EF_species = EFs[species_row[species_num]];
        writers = [];
        for writer_type in range(3):
            writers.append(setup_writer(data_types[writer_type], species_used[species_num], units[writer_type]));
        #calculate and write emissions for this species for each year 1997 - 2014
        for year in range(start_year, end_year+1):
            emissions_table = calculate_species_for_year(directory, species_num, EF_species, year, 0);
            regional_totals[year - start_year] += emissions_table;
            species_totals[year - start_year][0:7, species_num] = emissions_table[0:7, 14];
            plot_and_write_table(emissions_table, writers, plotter, species_num, EF_species, species_used[species_num] +"_" +  str(year));
        #do el nino years separately -- calculate and write emissions for July 1997 - June 1998
        el_nino_97_table = calculate_species_for_year(directory, species_num, EF_species, 1997, 7);
        plot_and_write_table(el_nino_97_table, writers, plotter, species_num, EF_species, species_used[species_num] + "_1997-1998 El Nino");
        la_nina_98_table = calculate_species_for_year(directory, species_num, EF_species, 1998, 7);        
        plot_and_write_table(la_nina_98_table, writers, plotter, species_num, EF_species, species_used[species_num] + "_1998-1999 La Nina");
        plot_and_write_table(el_nino_97_table - la_nina_98_table, writers, plotter, species_num, EF_species, species_used[species_num] + "_difference between 97-98 El Nino and 98-99 La Nina");
        print species_used[species_num] + " done";

    #calculate total emissions by adding up the results from each species, for each year
    print "plotting totals...";
    for regional_total in regional_totals:
        plot_regions_table(regional_total, plotter, species_num, EF_species, str(year) + " regional totals");
    for species_total in species_totals:
        plot_species_table(species_total, plotter, species_num, EF_species, str(year) + " all species");

    


# please compare this to http://www.falw.vu/~gwerf/GFED/GFED4/tables/GFED4.1s_CO.txt
if __name__ == '__main__':
        calculate_emissions()

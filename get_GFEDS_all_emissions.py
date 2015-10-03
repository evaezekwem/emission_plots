from plot import Plotter
import csv
import numpy as np
import h5py # if this creates an error please make sure you have the h5py library

def calculate_emissions():
    months       = '01','02','03','04','05','06','07','08','09','10','11','12'
    regions      = ['BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST', 'Global'];
    sources      = 'SAVA','BORF','TEMF','DEFO','PEAT','AGRI', 'All sources';
    species_used = 'CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3'
    species_row  = 2, 4, 13, 14, 3, 12, 8, 7, 33
    scar_values = [84, 4589, 270419, 41968, 632, 68299, 36987, 67141, 24906];

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

    start_year = 1997
    end_year   = 2014
    plotter = Plotter();

    for species_num in range(9):
        EF_species = EFs[species_row[species_num]];
        emissions_writer = csv.writer(open("./plots/tables/emissions/" + species_used[species_num] + '_emissions.csv', 'w'));
        scar_writer = csv.writer(open("./plots/tables/SCAR/" + species_used[species_num] + '_scar.csv','w'));
        scar_writer.writerow(["units - 2007 US$"]);
        emissions_writer.writerow(["units - 1E12 g"]);
        for year in range(start_year, end_year+1):
            emissions_table = np.zeros((7, 15)) # source, region
            source_emissions = np.zeros(7);
            emissions_writer.writerow([year]);
            emissions_writer.writerow([""] + regions);
            scar_writer.writerow([year]);
            scar_writer.writerow([""] + regions);
        
            string = directory+'/GFED4.1s_'+str(year)+'.hdf5'
            f = h5py.File(string, 'r')
            if year == start_year: # these are time invariable    
                basis_regions = f['/ancill/basis_regions'][:]
                grid_area     = f['/ancill/grid_cell_area'][:]
            
            total_emissions = np.zeros((720, 1440));
            for source in range(6):
                source_emissions = np.zeros((720, 1440))
                for month in range(12):
                    # read in DM emissions
                    string = '/emissions/'+months[month]+'/DM'
                    DM_emissions = f[string][:]
                    # read in the fractional contribution of each source
                    string = '/emissions/'+months[month]+'/partitioning/DM_'+sources[source]
                    contribution = f[string][:]
                    source_emissions += DM_emissions * contribution * EF_species[source];
                # calculate CO emissions as the product of DM emissions (kg DM per 
                    # m2 per month), the fraction the specific source contributes to 
                    # this (unitless), and the emission factor (g CO per kg DM burned)
                total_emissions += source_emissions;
                for region in range(15):
                    if region == 14:
                        mask = np.ones((720, 1440))
                    else:
                        mask = basis_regions == (region + 1)            
                    emissions_table[source, region] = np.sum(grid_area * mask * source_emissions);


            
            # fill table with total values for the globe (row 15) or basisregion (1-14)
            for region in range(15):
                if region == 14:
                    mask = np.ones((720, 1440))
                else:
                    mask = basis_regions == (region + 1)            
                emissions_table[6, region] = np.sum(grid_area * mask * total_emissions);
            print "\t" + str(year) + " done";
            # convert to $ value
            GRAMS_PER_TON = 907185;
            scar_table = emissions_table * scar_values[species_num] / GRAMS_PER_TON;
            # convert to Tg CO 
            emissions_table = emissions_table / 1E12;
            plotter.plot(species_used[species_num], year, "emissions", emissions_table);
            plotter.plot(species_used[species_num], year, "SCAR", scar_table);
            # write
            for source in range(7):
                emissions_source_list = emissions_table[source, :].tolist();
                scar_source_list = scar_table[source, :].tolist();
                emissions_source_list.insert(0, sources[source]);
                scar_source_list.insert(0, sources[source]);
                emissions_writer.writerow(emissions_source_list);
                scar_writer.writerow(scar_source_list);
                #writer.writerow([sources[source]] + emissions_table[source, :].tolist());
            emissions_writer.writerow([]);
            scar_writer.writerow([]);
            print emissions_table
        print species_used[species_num] + " done";

# please compare this to http://www.falw.vu/~gwerf/GFED/GFED4/tables/GFED4.1s_CO.txt
if __name__ == '__main__':
        calculate_emissions()

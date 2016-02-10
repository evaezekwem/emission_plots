from bokeh.charts import HeatMap, Bar, output_file, show
from bokeh.plotting import figure, output_file, save
from bokeh.palettes import Reds9 as red_palette
import numpy as np;
import pandas;

class Plotter:
    TOOLS="resize,crosshair,tap,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,hover";
    species      = ['CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3'];
    regions      = ['BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST'];
    sources      = ['SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI', 'All sources'];
    basic_sources      = ['SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI'];
    reindexed_sources = ['TEMF','BORF','SAVA','DEFO','PEAT','AGRI', 'All sources'];
    colors = ["#660033", "#FF0066", "#FFCC99", "#CCCC00", "#333300", "#00FF00", "#009999", "#66FFFF", "#000099", "#6600CC", "#CC7A00", "#FF9999", "#FFFF00", "#522900", "#006600"];
    NUM_YEARS = 18;
    NUM_SOURCES = 6;
    start_year = 1997;
    
    def format_table(self, table, column):
        formatted_data = pandas.DataFrame(data=table, index=self.sources, columns=column);
        formatted_data = formatted_data.reindex(self.reindexed_sources);
        return formatted_data;
    
    def format_time_series(self, time_series, column):
        time_series_sources = [];
        time_series_reindexed = [];
        for year in range(self.NUM_YEARS):
            for source in range(self.NUM_SOURCES + 1):
                time_series_sources.append(str(year+self.start_year) + " " + self.sources[source]);
                time_series_reindexed.append(str(year+self.start_year) + " " + self.reindexed_sources[source]);
        formatted_data = pandas.DataFrame(data=time_series, index=time_series_sources, columns=column);
        formatted_data = formatted_data.reindex(time_series_reindexed);
        return formatted_data;
        
    def format_source_time_series(self, time_series, column):
        time_series_sources = [];
        for year in range(self.NUM_YEARS):
            time_series_sources.append(str(year+self.start_year));
        formatted_data = pandas.DataFrame(data=time_series, index=time_series_sources, columns=column);
        return formatted_data;        
        
    def plot(self, identifier, chart, data, width):
        chart_title = chart + " - " + identifier;
        y_label = "Social Cost (2007 US $)";
        if(chart == "emissions"):
            y_label = "Total emissions (1E12 g)";
        output_file("plots/tables/" + chart + "/plots/" + identifier + "_" + chart + "_visualization.html", title = chart_title);
        
        p = Bar(data, width=width, height=700, title=chart_title, xlabel="Source", ylabel=y_label, legend="top_left", stacked=True, palette=self.colors, tools=self.TOOLS);
        save(p);
    
    def plot_heatmap(self, identifier, chart, average_table):
        chart_title = chart + " Heatmap";
        removed_total = np.delete(average_table, self.NUM_SOURCES, 0);
        formatted_data = None
        if(identifier == "regions"):
            formatted_data = pandas.DataFrame(data=removed_total, index=self.basic_sources, columns=self.regions);
        if(identifier == "species"):
            formatted_data = pandas.DataFrame(data=removed_total, index=self.basic_sources, columns=self.species);
        y_label = identifier;
        palette = red_palette[::-1]  # Reverse the color order so dark red is highest
        output_file("plots/tables/" + chart + "/plots/" + identifier + "_heatmap.html", title = chart_title)
        
        p = HeatMap(formatted_data, width=900, height=900, title=chart_title, xlabel="Source", ylabel=y_label, palette=palette);
	save(p);
   
    def plot_regions_total(self, identifier, chart, regions_table):
        formatted_data = self.format_table(regions_table, self.regions);
        self.plot(identifier, chart, formatted_data, 900);

    def plot_species_total(self, identifier, chart, species_table):
        formatted_data = self.format_table(species_table, self.species);
        self.plot(identifier, chart, formatted_data, 900);

    def plot_regions_time_series(self, identifier, chart, time_series):
        formatted_data = self.format_time_series(time_series, self.regions);
        self.plot(identifier, chart, formatted_data, 12000);
    
    def plot_species_time_series(self, identifier, chart, time_series):
        formatted_data = self.format_time_series(time_series, self.species);
        self.plot(identifier, chart, formatted_data, 12000);
        
    def plot_regions_source_time_series(self, identifier, chart, time_series):
        formatted_data = self.format_source_time_series(time_series, self.regions);
        self.plot(identifier, chart, formatted_data, 4000);
        


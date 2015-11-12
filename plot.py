from bokeh.charts import Bar, output_file, show
from bokeh.plotting import figure, output_file, save
import numpy as np;
import pandas;

class Plotter:
    TOOLS="resize,crosshair,tap,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,hover";
    species      = ['CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3'];
    regions      = ['BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST'];
    sources      = ['SAVA', 'BORF', 'TEMF', 'DEFO', 'PEAT', 'AGRI', 'All sources'];
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
        formatted_data = pandas.DataFrame(data=table, index=time_series_sources, columns=column);
        formatted_data = formatted_data.reindex(time_series_reindexed);
        return formatted_data;
        
    def plot(self, identifier, chart, data):
        chart_title = chart + " - " + identifier;
        y_label = "Social Cost (2007 US $)";
        if(chart == "emissions"):
            y_label = "Total emissions (1E12 g)";
        output_file("plots/tables/" + chart + "/plots/" + identifier + "_" + chart + "_visualization.html", title = chart_title);
        
        p = Bar(formatted_data, width=900, height=700, title=chart_title, xlabel="Source", ylabel=y_label, legend="top_left", stacked=True, palette=self.colors, tools=self.TOOLS);
        save(p);
    
    def plot_regions_total(self, identifier, chart, regions_table):
        formatted_data = self.format_table(self, regions_table, self.regions);
        plot(identifier, chart, formatted_data);

    def plot_species_total(self, identifier, chart, species_table):
        formatted_data = self.format_table(self, species_table, self.species);
        plot(identifier, chart, formatted_data);

    def plot_regions_time_series(self, identifier, chart, time_series):
        formatted_data = self.format_time_series(self, time_series, self.regions);
        plot(identifier, chart, formatted_data);
    
    def plot_species_time_series(self, identifier, chart, time_series):
        formatted_data = self.format_time_series(self, time_series, self.species);
        plot(identifier, chart, formatted_data);
        


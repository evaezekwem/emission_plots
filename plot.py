from bokeh.charts import Bar, output_file, show
from bokeh.plotting import figure, output_file, save
import numpy as np;
import pandas;

class Plotter:
    TOOLS="resize,crosshair,tap,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,hover";
    species      = ['CO2', 'CH4', 'BC', 'SO2', 'CO', 'OC', 'N2O', 'NOx', 'NH3'];
    regions      = ['BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST'];
    sources      = ['SAVA','BORF', 'TEMF','DEFO','PEAT','AGRI', 'All sources'];
    colors = ["#660033", "#FF0066", "#FFCC99", "#CCCC00", "#333300", "#00FF00", "#009999", "#66FFFF", "#000099", "#6600CC", "#CC7A00", "#FF9999", "#FFFF00", "#522900", "#006600"];
    
    def plot_regions_total(self, identifier, chart, species_table):
        without_global = species_table[0:7, 0:14];
        formatted_data = pandas.DataFrame(data=without_global, index=sources, columns=regions);
        formatted_data = formatted_data.reindex(['TEMF','BORF','SAVA','DEFO','PEAT','AGRI', 'All sources']);

        chart_title = chart + " - " + identifier;
        y_label = "Social Cost (2007 US $)";
        if(chart == "emissions"):
            y_label = "Total emissions (1E12 g)";
        output_file("plots/tables/" + chart + "/plots/" + identifier + "_" + chart + "_visualization.html", title = chart_title);
        
        p = Bar(formatted_data, width=900, height=700, title=chart_title, xlabel="Source", ylabel=y_label, legend="top_left", stacked=True, palette=region_colors, tools=TOOLS);
        save(p);

    def plot_species_total(self, identifier, chart, species_table):
        formatted_data = pandas.DataFrame(data=species_table, index=sources, columns=species);
        col_list = list(formatted_data);
        formatted_data = formatted_data.reindex(['TEMF','BORF','SAVA','DEFO','PEAT','AGRI', 'All sources']);

        chart_title = chart + " - " + identifier;
        y_label = "Social Cost (2007 US $)";
        if(chart == "emissions"):
            y_label = "Total emissions (1E12 g)";
        output_file("plots/tables/" + chart + "/plots/" + identifier + "_" + chart + "_visualization.html", title = chart_title);
        
        p = Bar(formatted_data, width=900, height=700, title=chart_title, xlabel="Source", ylabel=y_label, legend="top_left", stacked=True, palette=species_colors, tools=TOOLS);
        save(p);


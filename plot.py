from bokeh.charts import Bar, output_file, show
from bokeh.plotting import figure, output_file, save
import numpy as np;
import pandas;

class Plotter:
    
    def plot(self, species, identifier, chart, species_table):
        TOOLS="resize,crosshair,tap,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,hover";
        regions      = ['BONA', 'TENA', 'CEAM', 'NHSA', 'SHSA', 'EURO', 'MIDE', 'NHAF', 'SHAF', 'BOAS', 'TEAS', 'SEAS', 'EQAS', 'AUST'];
        sources      = 'TEMF','SAVA','BORF','DEFO','PEAT','AGRI', 'All sources';
        region_colors = ["#660033", "#FF0066", "#FFCC99", "#CCCC00", "#333300", "#00FF00", "#009999", "#66FFFF", "#000099", "#6600CC", "#CC7A00", "#FF9999", "#FFFF00", "#522900", "#006600"];
        
        without_global = species_table[0:7, 0:14];
        formatted_data = pandas.DataFrame(data=without_global, index=sources, columns=regions);
        col_list = list(formatted_data);
        col_list[0], col_list[1], col_list[2] = col_list[2], col_list[1], col_list[0];
        formatted_data.columns = col_list;

        chart_title = "" + species + " " + chart + " - " + identifier;
        y_label = "Social Cost (2007 US $)";
        if(chart == "emissions"):
            y_label = "Total emissions (1E12 g)";
        output_file("plots/tables/" + chart + "/plots/" + species + "_" + identifier + "_" + chart + "_visualization.html", title = chart_title);
        
        p = Bar(formatted_data, title=chart_title, xlabel="Source", ylabel=y_label, legend="top_left", stacked=True, palette=region_colors, tools=TOOLS);
        save(p);

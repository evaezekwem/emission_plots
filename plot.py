from get_GFED4s_CO_emissions import COEmissions;
from bokeh.plotting import figure, output_file, show
import numpy as np;

start_year = 1997
end_year = 2014

years = [1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014];

emissions = COEmissions()
CO_table = emissions.getCOData()

# output to static HTML file
output_file("plots/CO.html", title="Carbon Monoxide Emissions")

TOOLS="resize,crosshair,tap,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,hover";

# create a new plot with a title and axis labels
p = figure(title="Carbon Monoxide", x_axis_label='Year', y_axis_label='CO Emissions', tools=TOOLS);

colors = ["#660033", "#FF0066", "#FFCC99", "#CCCC00", "#333300", "#00FF00", "#009999", "#66FFFF", "#000099", "#6600CC", "#CC7A00", "#FF9999", "#FFFF00", "#522900", "#006600"];

# add a line renderer with legend and line thickness
for region in range(15):
    p.circle(years, CO_table[region, :], color=colors[region], size=10, alpha=0.5);
    p.line(years, CO_table[region, :], legend="region " + str(region), line_color=colors[region], line_width=2);

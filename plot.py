from get_GFED4s_CO_emissions import COEmissions;
from bokeh.plotting import figure, output_file, show

start_year = 1997
end_year = 2014

years = [1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014];

emissions = COEmissions()
CO_table = emissions.getCOData()

# output to static HTML file
output_file("CO.html", title="Carbon Monoxide Emissions")

# create a new plot with a title and axis labels
p = figure(title="Carbon Monoxide", x_axis_label='year', y_axis_label='y')

# add a line renderer with legend and line thickness
for region in range(15):
    p.line(years, CO_table[region, :], legend="Temp.", line_width=2)

# show the results
show(p)

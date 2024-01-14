import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import math
import csv
from cycler import cycler
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
plt.ioff()

plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 800
plt.rcParams['figure.figsize'] = 3,1.5
plt.rcParams['axes.prop_cycle'] = cycler(color=['blue','orange','brown','pink','gray','green','red','purple','#4B0082', '#2E8B57', '#7CFC00', '#C71585', '#808000', 
                                                '#228B22', '#008B8B', '#9ACD32', '#FF69B4', '#FF00FF', '#8B4513', '#708090', '#FFD700', '#006400', '#FF7F50', '#0000FF', 
                                                '#DC143C', '#1E90FF', '#FF0000', '#4682B4', '#556B2F', '#FFA07A', '#483D8B', '#BDB76B', '#FF4500', '#8A2BE2', '#CD5C5C', 
                                                '#8FBC8F', '#000080', '#FFFF00', '#8B0000', '#00FA9A', '#FF8C00'])


plt.rcParams['xtick.major.pad'] = 0.5
plt.rcParams['ytick.major.pad'] = 0.5
plt.rcParams['xtick.major.width'] = 0.2  # X-axis major tick width
plt.rcParams['ytick.major.width'] = 0.2
plt.rcParams['axes.titlepad'] = 1
plt.rcParams['axes.labelpad'] = 1
plt.rcParams['axes.linewidth'] = 0.24

#Sets editable parameters for all charts
legendsize = 6
axistitlesize = 7
legendtitlesize = 7
ticksize = 6
scalefactor = 0.5
titlesize = 12
axesspinecolour = 'gray'
allaxesspines = ['bottom','left','right','top']
line_width = 0.2

#Sets the axes ticksizes
plt.rcParams['xtick.labelsize'] = ticksize * scalefactor  # Change 12 to your desired font size
plt.rcParams['ytick.labelsize'] = ticksize * scalefactor

#Links to the model's outputted data
filepath1 = r"../Coupled water-sediment model/output_bed_elevation.csv"
filepath2 = r"../Coupled water-sediment model/output_depth.csv"
filepath3 = r"../Coupled water-sediment model/output_data.csv"
filepath4 = r"../Coupled water-sediment model/output_free_surface_elevation.csv"
filepath5 = r"../Coupled water-sediment model/output_bed_elevation_comparison.csv"
filepath6 = r"../Coupled water-sediment model/output_free_surface_comparison.csv"

lower_y_lim = 0
upper_y_lim = 90

#Gets relevant data from the parameters file
df_data = pd.read_csv(filepath3)

n_timesteps = int(df_data.loc[df_data['Parameter:'] == 'Number of time steps', 'Value:'].values[0])
n_steps = int(df_data.loc[df_data['Parameter:'] == 'Number of spatial steps', 'Value:'].values[0])
parameter_case = int(df_data.loc[df_data['Parameter:'] == 'Parameter case', 'Value:'].values[0])
x_tot = int(df_data.loc[df_data['Parameter:'] == 'Total length', 'Value:'].values[0])

upper_x_lim = x_tot
n_value = 20000

#Reads the two outputted files
df = pd.read_csv(filepath5,header=None)
df2 = pd.read_csv(filepath6,header=None)

#print(df2)
#print(df2.iloc[3])


#Creates the figure
fig1, ax1 = plt.subplots(1)

#Plots the data on the main graph
ax1.plot(df.iloc[0], df.iloc[1], label='Initial bed elevation $z_b$',color = 'darkgrey',linewidth = line_width, linestyle='--')
ax1.plot(df.iloc[0], df.iloc[2], label='Final bed elevation $z_b$', color = 'red',linewidth = line_width)
ax1.plot(df2.iloc[0], df2.iloc[1], label='Initial free surface elevation $η$',color = 'darkturquoise',linewidth = line_width * 0.75)
#ax1.plot(df2.iloc[0], df2.iloc[2], label='Maximum free surface elevation $η$',color = 'blue',linewidth = line_width * 0.75)
#ax1.plot(df2.iloc[0], df2.iloc[2], label='Final free surface elevation $η$', color = 'darkblue',linewidth = line_width)
ax1.plot(df2.iloc[0], df2.iloc[2], label='Final free surface elevation $η$', color = 'darkblue',linewidth = line_width * 0.75, linestyle = 'dashed')

#Creates the axis limits
ax1.set_xlim(0, upper_x_lim)
ax1.set_ylim(lower_y_lim, upper_y_lim)

#Creates the inset
#ax_inset = inset_axes(ax1, width="30%", height="30%", loc='upper right')
width, height = "33%", "33%"
ax_inset = inset_axes(ax1, width=width, height=height, borderpad=0.8, loc = 'lower left')


#Plots the data on the inset
ax_inset.plot(df2.iloc[0], df2.iloc[1], label = '_nolegend_',color = 'darkturquoise', linewidth = 0.75*line_width)
#ax_inset.plot(df2.iloc[0], df2.iloc[2], label='_nolegend_',color = 'blue',linewidth = 0.75*line_width)
ax_inset.plot(df2.iloc[0], df2.iloc[2], label = '_nolegend_', color = 'darkblue',linewidth = 0.75*line_width, linestyle = 'dashed')
ax_inset.plot(df.iloc[0], df.iloc[1], label = '_nolegend_', color='darkgrey',linewidth = line_width, linestyle='--')
ax_inset.plot(df.iloc[0], df.iloc[2], label = '_nolegend_', color='red', linewidth = line_width)
ax_inset.xaxis.set_major_locator(plt.MultipleLocator(5000))
ax_inset.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax_inset.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)

# Set limits for the zoomed-in portion
ax_inset.set_xlim(23000, 38000)
ax_inset.set_ylim(11, 38)

#Adds a grid to the main chart
ax1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)

#Modifies the box around the inset and marks it on the chart
mark_inset(ax1, ax_inset, loc1=2, loc2=4, fc="none", ec="0.5",linewidth = 0.2)

#Labels the axes
ax1.set_xlabel('Position $x$ (m)',fontsize=axistitlesize*scalefactor)
ax1.set_ylabel('Elevation $z$ (m)',fontsize=axistitlesize*scalefactor)

#Adds the legend and makes space for it
fig1.legend(fontsize = legendsize * scalefactor, loc='upper center', title='Legend:',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=3,bbox_to_anchor=(0.5, 1.0))#, bbox_to_anchor=(0.5, -0.03))
fig1.subplots_adjust(top=0.8)

#Saves the figure
save_file_path = r"../Results/Graphical results/Uniform_flow.png"
fig1.savefig(save_file_path, dpi=1000)

fig1.show()

# Your condition
difference_threshold = 0.045
row_index = 2
previous_row_index = 1
threshold_value = 28500

# Find the starting column index based on the threshold value in row 0
start_column = (df.iloc[0] > threshold_value).idxmax()

# Subset the DataFrame to columns starting from the selected column
subset_df = df.iloc[:, start_column:]

# Find the first column meeting the condition
column_index = (subset_df.iloc[row_index] - df.iloc[previous_row_index, start_column:] < difference_threshold).idxmax()

# Get the corresponding value from row 0
corresponding_value = df.at[0, column_index]

print('Gravel was transported a distance of:', 0.001* (corresponding_value - 28500), 'km')


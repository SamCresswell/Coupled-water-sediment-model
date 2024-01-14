import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import math
import csv
from cycler import cycler
import sys
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from colour import Color
plt.ioff()

plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 800
plt.rcParams['figure.figsize'] = 3,1.4
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

line_rate = 500

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
upper_y_lim = 6


output_data_filepath = r'../Coupled water-sediment model/output_data.csv'

#Gets relevant data from the parameters file
df_data = pd.read_csv(output_data_filepath)
print(df_data)

#Reads the number of timesteps, number of steps and final time
n_timesteps = int(df_data.loc[df_data['Parameter:'] == 'Number of time steps', 'Value:'].values[0])
n_steps = int(df_data.loc[df_data['Parameter:'] == 'Number of spatial steps', 'Value:'].values[0])
final_time = int(df_data.loc[df_data['Parameter:'] == 'Final time', 'Value:'].values[0])
parameter_case = int(df_data.loc[df_data['Parameter:'] == 'Parameter case', 'Value:'].values[0])
total_length = int(df_data.loc[df_data['Parameter:'] == 'Total length', 'Value:'].values[0])

display_case_1 = 0




upper_x_lim = total_length
n_value = 20000

colour1 = Color("red")
colours = list(colour1.range_to(Color("darkviolet"),n_timesteps))
matplotlib_colors = [c.rgb for c in colours] 

df3 = pd.read_csv(filepath2,header=None,usecols=range(n_steps+5))

#Creates the figure
fig, ax = plt.subplots(1)

for i in range(2,n_timesteps):
    
    if i % line_rate == 0:
        ax.plot(df3.iloc[0, 1:],df3.iloc[i,1:],color=matplotlib_colors[i],linewidth = 0.15)
        
    else:
        continue

ax.axhline(y=1.301, color='black', linestyle='--', linewidth=0.5, label='Analytical solution',zorder = 10)
ax.axhline(y=5, color='black', linestyle='dotted', linewidth=0.3, label='Initial flow depth $h_0$',zorder = 10)



#Creates the axis limits
ax.set_xlim(0, upper_x_lim)
ax.set_ylim(lower_y_lim, upper_y_lim)


#Adds a grid to the main chart
ax.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)

ax.minorticks_on()
ax.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)



#Labels the axes
ax.set_xlabel('Horizontal position $x$ (m)',fontsize=axistitlesize*scalefactor)
ax.set_ylabel('Flow depth $h$ (m)',fontsize=axistitlesize*scalefactor)

#Adds the legend and makes space for it
fig.legend(fontsize = legendsize * scalefactor, loc='upper center',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=3,bbox_to_anchor=(0.5, 1.0))#, bbox_to_anchor=(0.5, -0.03))

norm = Normalize(vmin=0, vmax=len(df3) - 1)
scalar_mappable = ScalarMappable(norm=norm, cmap=plt.cm.colors.ListedColormap(matplotlib_colors))  # You can choose a different colormap

# Add a colorbar
cbar = plt.colorbar(scalar_mappable, ax=ax)
cbar.set_label('Time $t$ (s)', fontsize=legendtitlesize * scalefactor)
cbar.ax.tick_params(labelsize= ticksize*scalefactor)
fig.tight_layout()

#Saves the figure
save_file_path = r"../Results/Graphical results/Uniform_flow_2.png"
fig.savefig(save_file_path, dpi=1000)

fig.show()




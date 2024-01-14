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
plt.rcParams['xtick.minor.width'] = 0.1  # X-axis major tick width
plt.rcParams['ytick.minor.width'] = 0.1
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
filepath0 = r"../Coupled water-sediment model/output_bed_elevation.csv"
filepath2 = r"../Coupled water-sediment model/output_depth.csv"
filepath3 = r"../Coupled water-sediment model/output_data.csv"
filepath4 = r"../Coupled water-sediment model/output_free_surface_elevation.csv"
filepath5 = r"../Coupled water-sediment model/output_bed_elevation_comparison.csv"
filepath6 = r"../Coupled water-sediment model/output_free_surface_comparison.csv"

filepath1 = r"C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 1, Sandbank\No_limiter_bed_elevation.csv"
filepath7 = r"C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 1, Sandbank\output_bed_elevation.csv"


lower_y_lim = 0
upper_y_lim = 1.15

#Gets relevant data from the parameters file
df_data = pd.read_csv(filepath3)

n_timesteps = int(df_data.loc[df_data['Parameter:'] == 'Number of time steps', 'Value:'].values[0])
n_steps = int(df_data.loc[df_data['Parameter:'] == 'Number of spatial steps', 'Value:'].values[0])
parameter_case = int(df_data.loc[df_data['Parameter:'] == 'Parameter case', 'Value:'].values[0])
x_tot = int(df_data.loc[df_data['Parameter:'] == 'Total length', 'Value:'].values[0])
final_time = int(df_data.loc[df_data['Parameter:'] == 'Final time', 'Value:'].values[0])


upper_x_lim = 1000
final_time = 40000

#Creates the analytical solution plot
start_value = 0
end_value = 1000
step_size = 0.5

x_positions = np.arange(0, 1001, 1).tolist()

#print(x_positions)

z_b_analytical = []

x_analytical = []


#print(x_positions)


for x in x_positions:
    
    #print(x_new)
    
    if (x > 500):
        
        z_b = 0
        
    elif (x < 200):
        
        z_b = 0
        
    else:
        
        z_b = math.sin((math.pi * (x - 300))/(200))
        
    xch = 10 * (1-z_b/10)**(-4)
    
    x_new = x + xch * 23800/2000
    
    x_analytical.append(x_new)
        
    #z_b_new = 1 + np.sin(math.pi*(1000- x + 20000/(2000)*10*(1-1/10 * np.sin((math.pi * (x - 300))/(200))^2))/(200))**2
    z_b_analytical.append(z_b)


#print(z_b_analytical)

#sys.exit()

#Reads the two outputted files
df = pd.read_csv(filepath1,header=None)
df2 = pd.read_csv(filepath5,header=None)



# Find the index of the row where the value in the first column matches 'specific_value'
index_value = df.iloc[:, 0][df.iloc[:, 0] == 23800].index

#index_value_2 = df.iloc[:, 0][df.iloc[:, 0] == 23800].index



#print(index_value)

#final_line = int(n_timesteps + 2)


#Creates the figure
fig1, ax1 = plt.subplots(1)


#Plots the data on the main graph
ax1.plot(df.iloc[0], df.iloc[1], label='Initial bed elevation $z_b$',color = 'darkgrey',linewidth = line_width, linestyle='solid')
ax1.plot(df.iloc[0], df.iloc[index_value[0], :], label='Non-slope-limited bed elevation $z_b$ at = 23800 s',color = 'black',linewidth = line_width, linestyle='solid')
#ax1.plot(df2.iloc[0], df2.iloc[index_value2[0], :], label='Bed elevation $z_b$ at {} seconds'.format(23800),color = 'darkgrey',linewidth = line_width, linestyle='solid')
ax1.plot(df2.iloc[0], df2.iloc[2], label='Slope-limited bed elevation $z_b$ at = 23800 s', color = 'blue',linewidth = line_width)
#ax1.plot(df2.iloc[0], df2.iloc[2], label='Final bed elevation $z_b$', color = 'red',linewidth = line_width)
ax1.plot(x_analytical, z_b_analytical, label='Analytical solution for t = 23800 s',color = 'red',linewidth = line_width, linestyle='dashed')
ax1.minorticks_on()
ax1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax1.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax1.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)

#Creates the axis limits
ax1.set_xlim(200, 800)
ax1.set_ylim(lower_y_lim, upper_y_lim)


#Adds a grid to the main chart
ax1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)

plt.axhline(y=1, color='grey', linestyle='--', linewidth=line_width)

#Labels the axes
ax1.set_xlabel('Horizontal position (m) $x$ (m)',fontsize=axistitlesize*scalefactor)
ax1.set_ylabel('Bed elevation $z_b$ (m)',fontsize=axistitlesize*scalefactor)

#Adds the legend and makes space for it
fig1.legend(fontsize = legendsize * scalefactor, loc='upper center', title='Legend:',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=3,bbox_to_anchor=(0.5, 1.0))#, bbox_to_anchor=(0.5, -0.03))
fig1.subplots_adjust(top=0.85)

#Saves the figure
save_file_path = r"../Results/Graphical results/Analytical_solution_hump.png"
fig1.savefig(save_file_path, dpi=1000)

fig1.show()





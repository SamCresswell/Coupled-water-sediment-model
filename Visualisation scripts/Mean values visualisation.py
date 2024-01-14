import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
from cycler import cycler
import sys
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

from colour import Color
plt.ioff()

#Set this to 0 for mean values and 1 for instantaneous values
display_case = 1

labelpadsize = 0.1
line_width = 0.25

plt.rcdefaults()
plt.rcParams["font.family"] = "Arial"
plt.rcParams['lines.linewidth'] = 0.1
plt.rcParams['axes.linewidth'] = 0.2
plt.rcParams['figure.dpi'] = 800
#plt.rcParams['figure.figsize'] = 3,6
all_colors = [
    '#8B0000', '#FF4500', '#FFA500', '#FFD700', '#808000', '#2E8B57',
    '#4682B4', '#000080', '#4B0082', '#8B008B', '#800080', '#A52A2A', '#696969', '#000000'
]


# Adding more randomization to the color list
#random.shuffle(all_colors)


plt.rcParams['axes.prop_cycle'] = cycler(color=all_colors)

plt.rcParams['xtick.major.pad'] = .07
plt.rcParams['ytick.major.pad'] = .07
plt.rcParams['xtick.major.width'] = 0.2  # X-axis major tick width
plt.rcParams['ytick.major.width'] = 0.2
plt.rcParams['axes.titlepad'] = 1
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['ytick.major.size'] = 1

#Sets editable parameters for all charts
legendsize = 10
axistitlesize = 10
ticksize = 8
scalefactor = 0.25
titlesize = 15
axesspinecolour = 'gray'
allaxesspines = ['bottom','left','right','top']

plt.rcParams['xtick.labelsize'] = ticksize*scalefactor
plt.rcParams['ytick.labelsize'] = ticksize*scalefactor

if display_case == 0:
    filepath = "C:/Users/sam_c/OneDrive - University of Glasgow/Invididual project/Model/Week 9 - Sediment model/Coupled water-sediment model/Coupled water-sediment model/output_mean_values.csv"
if display_case == 1:
    filepath = "C:/Users/sam_c/OneDrive - University of Glasgow/Invididual project/Model/Week 9 - Sediment model/Coupled water-sediment model/Coupled water-sediment model/output_instantaneous_values.csv"

output_data_filepath = "C:/Users/sam_c/OneDrive - University of Glasgow/Invididual project/Model/Week 9 - Sediment model/Coupled water-sediment model/Coupled water-sediment model/output_data.csv"

#Gets relevant data from the parameters file
df_data = pd.read_csv(output_data_filepath)
print(df_data)

#Reads the number of timesteps, number of steps and final time
n_timesteps = int(df_data.loc[df_data['Parameter:'] == 'Number of time steps', 'Value:'].values[0])
n_steps = int(df_data.loc[df_data['Parameter:'] == 'Number of spatial steps', 'Value:'].values[0])
final_time = int(df_data.loc[df_data['Parameter:'] == 'Final time', 'Value:'].values[0])
parameter_case = int(df_data.loc[df_data['Parameter:'] == 'Parameter case', 'Value:'].values[0])
total_length = int(df_data.loc[df_data['Parameter:'] == 'Total length', 'Value:'].values[0])





#Creates the figure and 3D axes
fig, axs = plt.subplots(4, 2,figsize=(3, 3))
ax1 = axs[0, 0]
ax2 = axs[0, 1]
ax3 = axs[1, 0]
ax4 = axs[1, 1]
ax5 = axs[2, 0]
ax6 = axs[2, 1]
ax7 = axs[3, 0]
ax8 = axs[3, 1]

for i in range(4):
    for j in range(2):

        #axs[i, j].set_yticks(list(ax1.get_xticks()) + [0])
        axs[i, j].grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
        axs[i, j].grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
        axs[i, j].set_xlabel('Specific flowrate $q$ (m$^3$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
        

df = pd.read_csv(filepath, header=None,usecols=range(11))#names = custom_headers)


custom_headers = ['t','mean_q','mean_u','mean_h','mean_conc','mean_q_b','mean_u_b','mean_En','mean_D','mean_Rn','mean_Courant_number']

df.columns = custom_headers

#print(df)

ax1.plot(df['t'], df['mean_q'],linewidth=line_width,color = all_colors[4])
ax1.set_xlabel('Time $t$ (s)', fontsize=axistitlesize*scalefactor)
ax1.set_ylabel('Specific flowrate $q$ (m$^3$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
ax1.set_xlim(left = 0)
#ax1.set_ylim(bottom = 0)

ax2.plot(df['mean_q'], df['mean_h'],linewidth=line_width,color = all_colors[9])
ax2.set_ylabel('Height $h$ (m)', fontsize=axistitlesize*scalefactor)
ax2.set_ylim(bottom = 0)
ax2.set_xlim(left = 0)

ax3.plot(df['mean_q'], df['mean_u'],linewidth=line_width,color = all_colors[1])
ax3.set_ylabel('Velocity $u$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
#ax3.set_ylim(bottom = 0)
ax3.set_xlim(left = 0)

ax4.plot(df['mean_q'], df['mean_conc'],linewidth=line_width,color = all_colors[3])
ax4.set_ylabel('Volumetric concentration $c$', fontsize=axistitlesize*scalefactor)
#ax4.set_ylim(bottom = 0)
ax4.set_xlim(left = 0)

ax5.plot(df['mean_q'], df['mean_q_b'],linewidth=line_width,color = all_colors[6])
ax5.set_ylabel('Bed flux $q_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
#ax5.set_ylim(bottom = 0)
#ax5.set_xlim(left = 0)

ax6.plot(df['mean_q'], df['mean_u_b'],linewidth=line_width,color = all_colors[10])
ax6.set_ylabel('Bed velocity $u_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
#ax6.set_ylim(bottom = 0)
ax6.set_xlim(left = 0)

ax7.plot(df['mean_q'], df['mean_En'],linewidth=line_width,color = all_colors[2])
ax7.set_ylabel('Entrainment flux $E$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
#ax7.set_ylim(bottom = 0)
ax7.set_xlim(left = 0)

ax8.plot(df['mean_q'], df['mean_D'],linewidth=line_width,color = all_colors[7])
ax8.set_ylabel('Deposition flux $D$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
#ax8.set_ylim(bottom = 0)
ax8.set_xlim(left = 0)

fig.tight_layout()
fig.subplots_adjust(hspace=.5)
plt.subplots_adjust(bottom=0.1)


ax1.set_title("a)", fontsize=titlesize * scalefactor, loc='left',y=1.03)
ax2.set_title("b)", fontsize=titlesize * scalefactor, loc='left',y=1.05)
ax3.set_title("c)", fontsize=titlesize * scalefactor, loc='left',y=1.05)
ax4.set_title("d)", fontsize=titlesize * scalefactor, loc='left',y=1.05)
ax5.set_title("e)", fontsize=titlesize * scalefactor, loc='left',y=1.05)
ax6.set_title("f)", fontsize=titlesize * scalefactor, loc='left',y=1.05)
ax7.set_title("g)", fontsize=titlesize * scalefactor, loc='left',y=1.05)
ax8.set_title("h)", fontsize=titlesize * scalefactor, loc='left',y=1.05)

save_file_path = r"../Results/Graphical results/Mean_values_single_slope.png"
plt.savefig(save_file_path, dpi=1000)

#fig.show()
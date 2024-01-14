import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
from cycler import cycler
import sys
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import random

from colour import Color
plt.ioff()

#Set this to 0 for mean values and 1 for instantaneous values
display_case = 3

labelpadsize = 0.1
line_width = 0.25

plt.rcParams["font.family"] = "Arial"
plt.rcParams['lines.linewidth'] = 0.1
plt.rcParams['axes.linewidth'] = 0.2
plt.rcParams['figure.dpi'] = 800
plt.rcParams['figure.figsize'] = 3*0.75,3*0.75
plt.rcParams['axes.titlepad'] = 0.4
plt.rcParams['axes.labelpad'] = 0.4


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

plt.rcParams['xtick.labelsize'] = ticksize * scalefactor  # Change 12 to your desired font size
plt.rcParams['ytick.labelsize'] = ticksize * scalefactor

cool_colors = ['#20B2AA', '#00CED1', '#48D1CC', '#87CEEB', '#1E90FF',
               '#7FFFD4', '#5F9EA0', '#40E0D0', '#87CEFA', '#4682B4',
               '#4169E1', '#6495ED', '#B0C4DE', '#AFEEEE', '#00BFFF',
               '#7B68EE', '#0000FF', '#87CEFA', '#AFEEEE', '#ADD8E6']

all_colors = [
    '#8B0000', '#FF4500', '#FFA500', '#FFD700', '#808000', '#2E8B57',
    '#4682B4', '#000080', '#4B0082', '#8B008B', '#800080', '#A52A2A', '#696969', '#000000'
]


# Adding more randomization to the color list
#random.shuffle(all_colors)


plt.rcParams['axes.prop_cycle'] = cycler(color=all_colors)

plt.rcParams['xtick.major.pad'] = 1
plt.rcParams['ytick.major.pad'] = 1
plt.rcParams['xtick.major.width'] = 0.2  # X-axis major tick width
plt.rcParams['ytick.major.width'] = 0.2
plt.rcParams['axes.titlepad'] = 1
plt.rcParams['axes.labelpad'] = 1
plt.rcParams['xtick.major.size'] = 1  # Length of major x-axis ticks
plt.rcParams['xtick.minor.size'] = 1  # Length of minor x-axis ticks
plt.rcParams['ytick.major.size'] = 1  # Length of major y-axis ticks
plt.rcParams['ytick.minor.size'] = 1  # Length of minor y-axis ticks



#Sets editable parameters for all charts
legendsize = 10
axistitlesize = 10
ticksize = 8
scalefactor = 0.25
titlesize = 15
legendtitlesize = 12
axesspinecolour = 'gray'
allaxesspines = ['bottom','left','right','top']

plt.rcParams['xtick.labelsize'] = ticksize*scalefactor
plt.rcParams['ytick.labelsize'] = ticksize*scalefactor

filepath = r"C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Results.xlsx"

output_data_filepath = "C:/Users/sam_c/OneDrive - University of Glasgow/Invididual project/Model/Week 9 - Sediment model/Coupled water-sediment model/Coupled water-sediment model/output_data.csv"



df = pd.read_excel(filepath, sheet_name='Thresholds',nrows=21)

#print(df)

fig1, axs = plt.subplots(2, 2)

ax1 = axs[0,0]
ax2 = axs[0,1]
ax3 = axs[1,0]
ax4 = axs[1,1]

ax1.plot(df['hyp'], df['u'],linewidth=line_width, color = 'blue')
ax2.plot(df['hyp'], df['q'],linewidth=line_width,color = 'red')
ax1.scatter(df['hyp'], df['u'], marker='x', color = 'blue', s = 3)
ax2.scatter(df['hyp'], df['q'], marker='x', color = 'red', s = 3)
ax3.plot(df['hyp'], df['u_shear'],linewidth=line_width,color = 'brown')
ax3.scatter(df['hyp'], df['u_shear'], marker='x', color = 'brown', s = 3)
ax4.plot(df['hyp'], df['tau_bx'],linewidth=line_width,color = 'green')
ax4.scatter(df['hyp'], df['tau_bx'], marker='x', color = 'green', s = 3)


ax2.set_xlabel('Volumetric concentration of \n simulated hyperconcentration (%)',fontsize=axistitlesize*scalefactor)
ax2.set_ylabel('Specific flowrate at initial bed mobilisation $q$ (ms$^-$$^2$)',fontsize=axistitlesize*scalefactor)

ax1.set_xlabel('Volumetric concentration of \n simulated hyperconcentration (%)',fontsize=axistitlesize*scalefactor)
ax1.set_ylabel('Depth-averaged velocity at initial bed mobilisation $u$ (ms$^-$$^1$)',fontsize=axistitlesize*scalefactor)

ax3.set_xlabel('Volumetric concentration of \n simulated hyperconcentration (%)',fontsize=axistitlesize*scalefactor)
ax3.set_ylabel('Shear velocity at initial bed mobilisation $u_*$ (ms$^-$$^1$)',fontsize=axistitlesize*scalefactor)

ax4.set_xlabel('Volumetric concentration of \n simulated hyperconcentration (%)',fontsize=axistitlesize*scalefactor)
ax4.set_ylabel('Bed stress at initial bed mobilisation $τ_{b{x}}$ (Nm$^-$$^2$)',fontsize=axistitlesize*scalefactor)

ax1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
ax2.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax2.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
ax3.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax3.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
ax4.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
ax4.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)

ax1.set_title("a)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
ax2.set_title("b)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
ax3.set_title("c)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
ax4.set_title("d)", fontsize=titlesize * scalefactor, loc='left',y=1.02)

#fig1.subplots_adjust(bottom=0.1)
fig1.tight_layout()

#fig1.show()

save_file_path = r"../Results/Graphical results/Mobilisation_thresholds.png"
fig1.savefig(save_file_path, dpi=1000)

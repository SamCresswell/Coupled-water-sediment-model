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
display_case = 1

labelpadsize = 0.1
line_width = 0.25

plt.rcParams["font.family"] = "Arial"
plt.rcParams['lines.linewidth'] = 0.1
plt.rcParams['axes.linewidth'] = 0.2
plt.rcParams['figure.dpi'] = 800
plt.rcParams['figure.figsize'] = 3,3.5

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

filepath = r"../Results/Processed results.xlsx"

output_data_filepath = r'../Coupled water-sediment model/output_data.csv'

"""
#Gets relevant data from the parameters file
df_data = pd.read_csv(output_data_filepath)
print(df_data)

#Reads the number of timesteps, number of steps and final time
n_timesteps = int(df_data.loc[df_data['Parameter:'] == 'Number of time steps', 'Value:'].values[0])
n_steps = int(df_data.loc[df_data['Parameter:'] == 'Number of spatial steps', 'Value:'].values[0])
final_time = int(df_data.loc[df_data['Parameter:'] == 'Final time', 'Value:'].values[0])
parameter_case = int(df_data.loc[df_data['Parameter:'] == 'Parameter case', 'Value:'].values[0])
total_length = int(df_data.loc[df_data['Parameter:'] == 'Total length', 'Value:'].values[0])
"""

df = pd.read_excel(filepath, sheet_name='Single_slope')

#print(df)

unique_values = df['case'].unique()

#print(unique_values)

mod_coeffs = df['coeff_3'].unique()

#print(df.columns)

#print(df['coeff_3'])

#print(mod_coeffs)

#sys.exit()

#Creates the figure and 3D axes
fig1, axs = plt.subplots(5, 2)
ax1 = axs[0, 0]
ax2 = axs[0, 1]
ax3 = axs[1, 0]
ax4 = axs[1, 1]
ax5 = axs[2, 0]
ax6 = axs[2, 1]
ax7 = axs[3, 0]
ax8 = axs[3, 1]
ax9 = axs[4, 0]

for i in range(5):
    for j in range(2):

        #axs[i, j].set_yticks(list(ax1.get_xticks()) + [0])
        axs[i, j].grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
        axs[i, j].grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
        axs[i, j].set_xlabel('Specific flowrate $q$ (m$^3$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
        
ax1.set_xlim(left = 0)
ax1.set_ylim(bottom = 0)

colours = ['red','green','blue','yellow','brown','orange','purple','grey','gold','magenta','lightblue','darkgrey','red','green','blue','yellow','brown','orange','purple','grey','gold','magenta','lightblue','darkgrey']

colourindex = 1

filtered_df = df[df['coeff_3'] == 1]

# Perform an action on each subset based on unique values
for value in unique_values:
    
    subset = df[df['case'] == value]

    
    #sys.exit()
    #current_colour = colours[colourindex]
    
    ax1.plot(subset['q'], subset['h'],linewidth=line_width)
    ax1.set_ylabel('Time $t$ (s)', fontsize=axistitlesize*scalefactor)
    ax1.set_ylabel('Height $h$ (m)', fontsize=axistitlesize*scalefactor)
    ax1.set_xlim(left = 0)
    #ax1.set_ylim(bottom = 0)

    ax2.plot(subset['q'], subset['u'],linewidth=line_width)#,color = current_colour)
    ax2.set_ylabel('Velocity $u$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    #ax2.set_ylim(bottom = 0)
    ax2.set_xlim(left = 0)
    
    ax3.plot(subset['q'], subset['conc'],linewidth=line_width)#,color = current_colour)
    ax3.set_ylabel('Volumetric concentration $c$', fontsize=axistitlesize*scalefactor)
    #ax3.set_ylim(bottom = 0)
    ax3.set_xlim(left = 0)

    ax4.plot(subset['q'], subset['q_b'],linewidth=line_width)#,color = current_colour)
    ax4.set_ylabel('Bedload flux $q_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    #ax4.set_ylim(bottom = 0)
    ax4.set_xlim(left = 0)

    ax5.plot(subset['q'], subset['u_b'],linewidth=line_width)#,color = current_colour)
    ax5.set_ylabel('Bed velocity $u_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    #ax5.set_ylim(bottom = 0)
    ax5.set_xlim(left = 0)

    ax6.plot(subset['q'], subset['E'],linewidth=line_width)#,color = current_colour)
    ax6.set_ylabel('Entrainment flux $E$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    #ax6.set_ylim(bottom = 0)
    ax6.set_xlim(left = 0)

    ax7.plot(subset['q'], subset['D'],linewidth=line_width)#,color = current_colour)
    ax7.set_ylabel('Deposition flux $D$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax7.set_ylim(bottom = 0)
    ax7.set_xlim(left = 0)

    ax8.plot(subset['q'], subset['R_n'],linewidth=line_width)#,color = current_colour)
    ax8.set_ylabel('Rouse number $R_n$', fontsize=axistitlesize*scalefactor)
    #ax8.set_ylim(bottom = 0)
    ax8.set_xlim(left = 0)

    ax9.plot(subset['q'], subset['Co'],linewidth=line_width)#,color = current_colour)
    ax9.set_ylabel('Courant number $C$', fontsize=axistitlesize*scalefactor)
    #ax9.set_ylim(bottom = 0)
    ax9.set_xlim(left = 0)    
    
    colourindex += 1



fig1.tight_layout()
#fig1.show()

#Creates the second figure and 3D axes
fig1b, axs = plt.subplots(2, 2,figsize=(2.5, 2.5))
axb1 = axs[0, 0]
axb2 = axs[0, 1]
axb3 = axs[1, 0]
axb4 = axs[1, 1]


for i in range(2):
    for j in range(2):

        #axs[i, j].set_yticks(list(ax1.get_xticks()) + [0])
        #axs[i, j].minorticks_on()
        axs[i, j].grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
        axs[i, j].grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
        #axs[i, j].grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
        axs[i, j].grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
        axs[i, j].set_xlabel('Specific flowrate $q$ (m$^3$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
        
ax1.set_xlim(left = 0)
ax1.set_ylim(bottom = 0)

colours = ['red','green','blue','yellow','brown','orange','purple','grey','gold','magenta','lightblue','darkgrey','red','green','blue','yellow','brown','orange','purple','grey','gold','magenta','lightblue','darkgrey']

colourindex = 1

filtered_df = df[df['coeff_3'] == 1]

unique_values2 = filtered_df['Hyp'].unique()

unique_coeffs = df['coeff_3'].unique().tolist()

# Perform an action on each subset based on unique values
for value in unique_values2:
    
    subset = filtered_df[filtered_df['Hyp'] == value]
    subset.reset_index(drop=True, inplace=True)

    axb1.plot(subset['q'], subset['conc'],linewidth=line_width,label = value)
    axb1.set_ylabel('Volumetric concentration $c$', fontsize=axistitlesize*scalefactor)
    axb1.set_xlim(left = 0, right = 80)
    axb1.set_yscale('log')
    axb1.tick_params(axis='y', width=line_width)
    axb1.tick_params(axis='y', which='minor', width=line_width)
    #axb1.set_ylim(bottom = 0.1,top=0.05)


    axb2.plot(subset['q'], subset['q_b'],linewidth=line_width, label = '_nolegend_')
    axb2.set_ylabel('Bed flux $q_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    axb2.set_ylim(bottom = 0,top = 0.06)
    axb2.set_xlim(left = 0, right = 80)
    
    axb3.plot(subset['q'], subset['u_b'],linewidth=line_width, label = '_nolegend_')
    axb3.set_ylabel('Bed velocity $u_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    axb3.set_ylim(bottom = 0,top = 0.26)
    axb3.set_xlim(left = 0, right = 80)

    axb4.plot(subset['q'], subset['E'],linewidth=line_width, label = '_nolegend_')
    axb4.set_ylabel('Entrainment flux $E$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    axb4.set_ylim(bottom = 0,top=0.007)
    axb4.set_xlim(left = 0, right = 80)


# Displaying only the first four items in the legend
fig1b.legend(fontsize = legendsize * scalefactor, loc='lower center', title='Volumetric concentration of simulated hyperconcentration (%):',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=6)
axb1.set_title("a)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
axb2.set_title("b)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
axb3.set_title("c)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
axb4.set_title("d)", fontsize=titlesize * scalefactor, loc='left',y=1.02)


fig1b.tight_layout()
fig1b.subplots_adjust(bottom=0.1)
#fig1b.show()
save_file_path = r"../Results/Graphical results/HYperconcentration_impacts.png"
fig1b.savefig(save_file_path, dpi=1000)





fig2 = plt.figure()
ax10 = fig2.add_subplot(1, 1, 1)

for value in unique_values:
    
    subset = df[df['case'] == value]

    if subset['coeff_3'].iloc[0] == 1:
    
        ax10.plot(subset['q'], subset['q_b'],marker='x',markersize = 0.2,linewidth=0.1,color = 'blue')
     
    if subset['coeff_3'].iloc[0] == 2:
    
        ax10.plot(subset['q'], subset['q_b'],marker='x',markersize = 0.2,linewidth=0.1,color = 'red')
        
    if subset['coeff_3'].iloc[0] == 4:
    
        ax10.plot(subset['q'], subset['q_b'],marker='x',markersize = 0.2,linewidth=0.1,color = 'green')
        
    if subset['coeff_3'].iloc[0] == 6:
    
        ax10.plot(subset['q'], subset['q_b'],marker='x',markersize = 0.2,linewidth=0.1,color = 'pink')
        
    ax10.set_ylabel('Bed flux $q_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax10.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
    ax10.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
    #ax10.set_xscale('log')
    
ax10.set_xlim(left = 0)
ax10.set_ylim(bottom = 0)
fig2.tight_layout()
#fig2.show()

fig2b, (ax2b1, ax2b2) = plt.subplots(1,2,figsize=(2.5, 1.5))


for value in (unique_coeffs):
    
    index = unique_coeffs.index(value)
    colour = all_colors[index]
    subset = df[df['coeff_3'] == value]
    subset1 = subset[subset['Hyp'] == 40]
    
    ax2b1.plot(subset1['q'], subset1['q_b'],marker='x',markersize = 0.2,linewidth=0.2,color = colour, label = value)
    ax2b2.plot(subset1['q'], subset1['u_b'],marker='x',markersize = 0.2,linewidth=0.2,color = colour, label = '_nolegend_',)
    
    
    subset2 = subset[subset['Hyp'] == 0]
    
    ax2b2.plot(subset2['q'], subset2['u_b'],marker='x',markersize = 0.2,linewidth=0.2,color = colour, label = '_nolegend_',linestyle = 'dotted')
    ax2b1.plot(subset2['q'], subset2['q_b'],marker='x',markersize = 0.2,linewidth=0.2,color = colour, label = '_nolegend_',linestyle = 'dotted')
    
    ax2b2.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
    ax2b2.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
    ax2b1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
    ax2b1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
    
    ax2b1.set_xlabel('Specific flowrate $q$ (m$^3$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax2b2.set_xlabel('Specific flowrate $q$ (m$^3$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax2b1.set_ylabel('Bed flux $q_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax2b2.set_ylabel('Bed velocity $u_b$ (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    
ax2b1.set_xlim(left = 0)
ax2b2.set_xlim(left = 0)
ax2b1.set_ylim(bottom = 0)
ax2b2.set_ylim(bottom = 0)

fig2b.tight_layout()
fig2b.subplots_adjust(bottom=0.2)

fig2b.legend(fontsize = legendsize * scalefactor, loc='lower center', title='Bedload flux modification coefficient $𝜑$:',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=4)
save_file_path = r"../Results/Graphical results/modification_coefficient.png"
fig2b.savefig(save_file_path, dpi=1000)


fig2b.show()


fig3 = plt.figure()
ax11 = fig3.add_subplot(1, 1, 1)

for value in unique_values:
    
    subset = df[df['case'] == value]

    if subset['coeff_3'].iloc[0] == 1:
    
        ax11.plot(subset['q'], subset['E'],marker='x',markersize = 0.2,linewidth=0.1,color = 'blue')
     
    if subset['coeff_3'].iloc[0] == 2:
    
        ax11.plot(subset['q'], subset['E'],marker='x',markersize = 0.2,linewidth=0.1,color = 'red')
        
    if subset['coeff_3'].iloc[0] == 4:
    
        ax11.plot(subset['q'], subset['E'],marker='x',markersize = 0.2,linewidth=0.1,color = 'green')
        
    ax10.set_ylabel('Entrainment flux $E$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax10.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
    ax10.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
    #ax10.set_xscale('log')
    
ax10.set_xlim(left = 0)
ax10.set_ylim(bottom = 0)
#fig.tight_layout()
#fig3.show()

fig4 = plt.figure()
ax12 = fig4.add_subplot(1, 1, 1)


for value in unique_values:
    
    subset = df[df['case'] == value]
    
    #print(subset.columns)
    #sys.exit()

    if subset['coeff_3'].iloc[0] == 1:
    
        ax12.plot(subset['q'], subset['conc'],marker='x',markersize = 0.2,linewidth=0.1,color = 'blue')
     
    if subset['coeff_3'].iloc[0] == 2:
    
        ax12.plot(subset['q'], subset['conc'],marker='x',markersize = 0.2,linewidth=0.1,color = 'red')
        
    if subset['coeff_3'].iloc[0] == 4:
    
        ax12.plot(subset['q'], subset['conc'],marker='x',markersize = 0.2,linewidth=0.1,color = 'green')
        
    ax12.set_ylabel('Volumetric concentration $c$', fontsize=axistitlesize*scalefactor)
    ax12.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
    ax12.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
    #ax10.set_xscale('log')
    
ax12.set_xlim(left = 0)
ax12.set_ylim(bottom = 0)
fig4.tight_layout()
save_file_path = r'C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Diagrams\Graphical outputs\Case 6, Single slope\Concentration_graph.png'  # Replace this with your desired file path
#fig4.savefig(save_file_path, dpi=1000)
#fig4.show()

fig5 = plt.figure()
ax13 = fig3.add_subplot(1, 1, 1)

for value in unique_values:
    
    subset = df[df['case'] == value]

    if subset['coeff_3'].iloc[0] == 1:
    
        ax13.plot(subset['q'], subset['E'],marker='x',markersize = 0.2,linewidth=0.1)
     
        
    ax13.set_ylabel('Entrainment flux $E$ (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)
    ax13.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
    ax13.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)
    #ax10.set_xscale('log')
    
ax13.set_xlim(left = 0)
ax13.set_ylim(bottom = 0)
fig5.tight_layout()




#fig5.show()

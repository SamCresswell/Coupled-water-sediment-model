import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler
import sys
plt.ioff()

plt.rcdefaults()
plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 800
plt.rcParams['figure.figsize'] = 3,3
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

all_colours = [
    '#8B0000', '#FF4500', '#FFA500', '#FFD700', '#808000', '#2E8B57',
    '#4682B4', '#000080', '#4B0082', '#8B008B', '#800080', '#A52A2A', '#696969', '#000000'
]

plt.rcParams['xtick.labelsize'] = ticksize * scalefactor  # Change 12 to your desired font size
plt.rcParams['ytick.labelsize'] = ticksize * scalefactor

filepath = r"../Results/Processed results.xlsx"

# Read the Excel file into a pandas DataFrame
df = pd.read_excel(filepath, sheet_name='Particle_diameter')

print(df)
#sys.exit()

fig1, axs = plt.subplots(2,2)

ax1 = axs[0,0]
ax2 = axs[0,1]
ax3 = axs[1,0]
ax4 = axs[1,1]

#filtered_df = df[df['Hyp'] != 40]

unique_values= df['d_median'].unique().tolist()

for value in unique_values:
    
    index = unique_values.index(value)
    
    colour = all_colours[index]
    
    subset = df[df['d_median'] == value]
    
    ax1.plot(subset['q'], subset['u'],linewidth=line_width, color = colour,zorder = 10, markersize=1, marker = 'x', markeredgewidth = .1, label = '_nolegend_')
    ax2.plot(subset['q'], subset['conc'],linewidth=line_width, color = colour,zorder = 10, markersize=1, marker = 'x', markeredgewidth = .1, label = '_nolegend_')
    ax3.plot(subset['q'], subset['q_b'],linewidth=line_width, color = colour,zorder = 10, markersize=1, marker = 'x', markeredgewidth = .1, label = value * 1000)
    ax4.plot(subset['q'], subset['u_b'],linewidth=line_width, color = colour,zorder = 10, markersize=1, marker = 'x', markeredgewidth = .1, label = '_nolegend_')

ax1.set_xlabel('Specific flowrate $q$ (m$^2$s$^-$$^1$)',fontsize=axistitlesize*scalefactor)
ax2.set_xlabel('Specific flowrate $q$ (m$^2$s$^-$$^1$)',fontsize=axistitlesize*scalefactor)
ax3.set_xlabel('Specific flowrate $q$ (m$^2$s$^-$$^1$)',fontsize=axistitlesize*scalefactor)
ax4.set_xlabel('Specific flowrate $q$ (m$^2$s$^-$$^1$)',fontsize=axistitlesize*scalefactor)

ax1.set_title("a)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
ax2.set_title("b)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
ax3.set_title("c)", fontsize=titlesize * scalefactor, loc='left',y=1.02)
ax4.set_title("d)", fontsize=titlesize * scalefactor, loc='left',y=1.02)

ax1.set_ylabel('Flow velocity $u$ (ms$^{-1}$)',fontsize=axistitlesize*scalefactor)
ax2.set_ylabel('Gravel volumetric concentration in flow $c$',fontsize=axistitlesize*scalefactor)
ax3.set_ylabel('Bedload flux $q_b$ (m$^2$s$^-$$^1$)',fontsize=axistitlesize*scalefactor)
ax4.set_ylabel('Bed velocity $u_b$ (ms$^-$$^1$)',fontsize=axistitlesize*scalefactor)

#ax1.set_ylim(bottom = 0)
ax2.set_yscale('log')
ax1.tick_params(axis='y', width=line_width)
ax1.tick_params(axis='y', which='minor', width=line_width)

ax1.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
ax4.minorticks_on()

ax1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax1.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax1.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax2.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax2.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax2.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax2.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax3.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax3.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax3.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax3.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax4.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax4.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax4.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax4.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax1.tick_params(which='minor', width=0.1)
ax2.tick_params(which='minor', width=0.1)
ax3.tick_params(which='minor', width=0.1)
ax4.tick_params(which='minor', width=0.1)
ax1.tick_params(which='minor', width=0.1)
ax2.set_ylim(bottom = 0)
ax3.set_ylim(bottom = 0)
ax4.set_ylim(bottom = 0)
ax1.set_xlim(left = 0)
ax2.set_xlim(left = 0)
ax3.set_xlim(left = 0)
ax4.set_xlim(left = 0)

fig1.legend(fontsize = legendsize * scalefactor, loc='upper center', title='Bed median grain size $d_{50}$ (mm):',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=6)

#fig1.show()
fig1.tight_layout()
fig1.subplots_adjust(top=0.92)

save_file_path = r"../Results/Graphical results/Particle_grain_size_effects.png"
fig1.savefig(save_file_path, dpi=1000)

import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler
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

all_colours = [
    '#8B0000', '#FF4500', '#FFA500', '#FFD700', '#808000', '#2E8B57',
    '#4682B4', '#000080', '#4B0082', '#8B008B', '#800080', '#A52A2A', '#696969', '#000000'
]

plt.rcParams['xtick.labelsize'] = ticksize * scalefactor  # Change 12 to your desired font size
plt.rcParams['ytick.labelsize'] = ticksize * scalefactor

filepath = r"../Results/Processed results.xlsx"

# Read the Excel file into a pandas DataFrame
df = pd.read_excel(filepath, sheet_name='Double_slope_entrainment',nrows = 45)

#print(df)

fig1, ax1 = plt.subplots(1)

unique_values = df['hyp'].unique().tolist()




for value in unique_values:
    
    index = unique_values.index(value)
    
    colour = all_colours[index]
    
    subset = df[df['hyp'] == value]
    
    ax1.plot(subset['q_peak'], subset['distance'],linewidth=line_width, color = colour,zorder = 10, markersize=5, marker = 'x', markeredgewidth = .1, label = value)
    

ax1.set_xlabel('Peak flowrate $q$ (m$^2$s$^-$$^1$)',fontsize=axistitlesize*scalefactor)
ax1.set_ylabel('Gravel deposition distance downstream from GST (km)',fontsize=axistitlesize*scalefactor)
ax1.set_ylim(bottom = 0)
ax1.set_xlim(left = 0)

ax1.minorticks_on()
ax1.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax1.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax1.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax1.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
plt.tick_params(which='minor', width=0.1)

fig1.legend(fontsize = legendsize * scalefactor, loc='upper center', title='Volumetric concentration of simulated hyperconcentration (%):',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=6)
fig1.subplots_adjust(top=0.85)

#fig1.show()

save_file_path = r"../Results/Graphical results/entrainment_distance.png"
fig1.savefig(save_file_path, dpi=1000)

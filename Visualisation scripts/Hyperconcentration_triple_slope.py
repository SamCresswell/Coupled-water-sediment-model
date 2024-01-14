import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler
import sys
import numpy as np
import matplotlib.patches as patches

#Set to 1 for no hyperconcentration, 2 for 4%, 3 for 8%, 4 for 12%
case_to_display = 1

q_peak = 30

    
hyp_conc = 0
save_file_path = r"../Results/Graphical results/Three_slope_model_30.png"
ups_x_lower_lim, ups_x_upper_lim = 0,1000
gst_x_lower_lim, gst_x_upper_lim = 24000, 38000
ds_x_lower_lim, ds_x_upper_lim = 45000, 50000
ups_y_lower_lim, ups_y_upper_lim = 76.5, 79
gst_y_lower_lim, gst_y_upper_lim = 12, 29
ds_y_lower_lim, ds_y_upper_lim = 0,6
text_pos_ups_x, text_pos_ups_y = 5000, 78
text_pos_gst_x, text_pos_gst_y = 30000, 8
text_pos_ds_x, text_pos_ds_y = 46000, 12
lower_x_lim = 0
upper_x_lim = 50000
lower_y_lim = -11
upper_y_lim = 82
   


plt.ioff()
plt.rcdefaults()
plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 800
#plt.rcParams['figure.figsize'] = 3,7
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
plt.rcParams['ytick.minor.width'] = 0.1
plt.rcParams['xtick.minor.width'] = 0.1

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

filepath1 = r"../Results/Final_bed_elevations.csv"
filepath2 = r"../Results/Final_free_surfaces.csv"
filepath3 = r"../Coupled water-sediment model/output_bed_elevation_comparison.csv"

# Read the Excel file into a pandas DataFrame
df1 = pd.read_csv(filepath1, header = None, usecols=range(207))
df2 = pd.read_csv(filepath1, header = None, usecols=range(207))
df3 = pd.read_csv(filepath3,header=None)

#Extracting the fourth to last column index
fourth_to_last_col_index = 4

first_row_values = df1.iloc[0, fourth_to_last_col_index:]

# Extracting the fourth to last column values of all rows (excluding the first row)
other_rows_values = df1.iloc[1:, fourth_to_last_col_index:]



#Creates several colourmapping sets
colormapping1 = {
    60: '#0000FF',  # Blue
    50: '#008000',  # Green
    40: '#ADD8E6',  # Light Blue
    30: '#90EE90',  # Light Green
    25: '#00008B',  # Dark Blue
    20: '#006400',  # Dark Green
    15: '#00FFFF',  # Cyan
}

# Colormap 2
colormapping2 = {
    60: '#800080',  # Purple
    50: '#FF0000',  # Red
    40: '#FFA500',  # Orange
    30: '#DAA520',  # Goldenrod
    25: '#4B0082',  # Indigo
    20: '#00FA9A',  # Medium Spring Green
    15: '#FF6347',  # Tomato
}

# Colormap 3
colormapping3 = {
    60: '#008000',  # Green
    50: '#0000FF',  # Blue
    40: '#FFFF00',  # Yellow
    30: '#800080',  # Purple
    25: '#FFA07A',  # Light Salmon
    20: '#BC8F8F',  # Rosy Brown
    15: '#00FF7F',  # Spring Green
}

# Colormap 4
colormapping4 = {
    60: '#FF0000',  # Red
    50: '#800080',  # Purple
    40: '#FFB6C1',  # Light Pink
    30: '#008080',  # Teal
    25: '#008B8B',  # Dark Cyan
    20: '#E9967A',  # Dark Salmon
    15: '#DDA0DD',  # Plum
}

colormapping5 = {
    0: '#8B0000',  # Dark Red
    4: '#FF4500',  # Orange Red
    8: '#FFA500',  # Orange
    12: '#FFD700',  # Gold
    25: '#808000',  # Olive
    20: '#2E8B57',  # Sea Green
    15: '#4682B4',  # Steel Blue
}

#Creating the figure and subplots
fig, axs = plt.subplots(3, 2,figsize = (2.5,3.2))

#Creates larger subplots spanning all three columns
big_subplot = plt.subplot2grid((3, 1), (0, 0), colspan=3)
ax2 = plt.subplot2grid((3, 1), (1, 0), colspan=3)

#Adds lines marking the initial position of the GST
big_subplot.axvline(x=26000, color='blue', linestyle='dashed', label='Initial start position of GST',linewidth = line_width)
big_subplot.axvline(x=31000, color='red', linestyle='dashed', label='Initial end position of GST',linewidth = line_width)

# Creating the big subplot in the first row for a non-hyperconcentrated flood
for index, row in other_rows_values.iterrows():
    #big_subplot.plot(first_row_values, row, label = df1.iloc[index, 0],linewidth = line_width)

    if int(df1.iloc[index, 0]) == q_peak:

        color = colormapping5.get(int(df1.iloc[index, 1]))
        
        big_subplot.plot(first_row_values, row, label = df1.iloc[index, 1]+' %',linewidth = line_width, color = color)
        
    big_subplot.minorticks_on()
    big_subplot.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
    big_subplot.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
    big_subplot.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
    big_subplot.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
        

#Adds a line showing the initial bed profile
big_subplot.plot(df3.iloc[0], df3.iloc[1], label='Initial',color = 'black',linewidth = line_width, linestyle='--')

#Defines the smaller plots showing the upstream and downstream ends of the model
ax1 = axs[2, 0]
ax3 = axs[2, 1]

# Creating insets in the second row
for i in range(2):
    
    axs[2, i].plot(df3.iloc[0], df3.iloc[1], label='_nolegend_',color = 'black',linewidth = line_width, linestyle='--',zorder= 10)

    
    for index, row in other_rows_values.iterrows():
        
        if int(df1.iloc[index, 0]) == q_peak:
    
            color = colormapping5.get(int(df1.iloc[index, 1]))
            
            axs[2, i].plot(first_row_values, row, label = '_nolegend_',linewidth = line_width, color = color)
            
    axs[2, i].minorticks_on()
    axs[2, i].grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
    axs[2, i].grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
    axs[2, i].grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
    axs[2, i].grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)

ax2.plot(df3.iloc[0], df3.iloc[1], label='_nolegend_',color = 'black',linewidth = line_width, linestyle='--',zorder= 10)
for index, row in other_rows_values.iterrows():
    
    if int(df1.iloc[index, 0]) == q_peak:

        color = colormapping5.get(int(df1.iloc[index, 1]))
        
        ax2.plot(first_row_values, row, label = '_nolegend_',linewidth = line_width, color = color)
        
ax2.minorticks_on()
ax2.grid(visible=True, axis='y', color='lightgray', linewidth=0.1, zorder = 0.1)
ax2.grid(visible=True, axis='x', color='lightgray', linewidth=0.1, zorder = 0.1)
ax2.grid(visible=True, axis='x', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
ax2.grid(visible=True, axis='y', which='minor', color='lightgray', linewidth=0.05, zorder = 0.1)
    
ax1.set_title('Upstream end',fontsize=titlesize*0.5*scalefactor,y=1)
ax2.set_title('Gravel-sand transition',fontsize=titlesize*0.5*scalefactor,y=1)
ax3.set_title('Downstream end',fontsize=titlesize*0.5*scalefactor,y=1)

ax1.set_xlim(ups_x_lower_lim, ups_x_upper_lim)
ax2.set_xlim(gst_x_lower_lim, gst_x_upper_lim )
ax3.set_xlim(ds_x_lower_lim, ds_x_upper_lim )

ax1.set_ylim(ups_y_lower_lim, ups_y_upper_lim )
ax2.set_ylim(gst_y_lower_lim, gst_y_upper_lim )
ax3.set_ylim(ds_y_lower_lim, ds_y_upper_lim )

rect1 = patches.Rectangle((ups_x_lower_lim , ups_y_lower_lim,), ups_x_upper_lim - ups_x_lower_lim, ups_y_upper_lim - ups_y_lower_lim, linewidth=line_width, edgecolor='grey', facecolor='none')
rect2 = patches.Rectangle((gst_x_lower_lim , gst_y_lower_lim,), gst_x_upper_lim - gst_x_lower_lim,gst_y_upper_lim - gst_y_lower_lim, linewidth=line_width, edgecolor='grey', facecolor='none')
rect3 = patches.Rectangle((ds_x_lower_lim , ds_y_lower_lim), ds_x_upper_lim - ds_x_lower_lim , ds_y_upper_lim - ds_y_lower_lim, linewidth=line_width, edgecolor='grey', facecolor='none')

big_subplot.add_patch(rect1)
big_subplot.add_patch(rect2)
big_subplot.add_patch(rect3)

ax2.axvline(x=26000, color='blue', linestyle='dashed', label='_nolegend_',linewidth = line_width)
ax2.axvline(x=31000, color='red', linestyle='dashed', label='_nolegend_',linewidth = line_width)

big_subplot.text(text_pos_ups_x, text_pos_ups_y , 'Upstream end', ha='center', va='center', color='black',fontsize=ticksize*scalefactor )
big_subplot.text(text_pos_gst_x, text_pos_gst_y , 'Gravel-sand transition', ha='center', va='center', color='black',fontsize=ticksize*scalefactor )
big_subplot.text(text_pos_ds_x, text_pos_ds_y , 'Downstream \n end', ha='center', va='center', color='black',fontsize=ticksize*scalefactor )

#Creates the axis limits
big_subplot.set_xlim(lower_x_lim, upper_x_lim)
big_subplot.set_ylim(lower_y_lim, upper_y_lim)

#Labels the axes
big_subplot.set_xlabel('Horizontal position $x$ (m)',fontsize=axistitlesize*scalefactor)
big_subplot.set_ylabel('Bed elevation $z_b$ (m)',fontsize=axistitlesize*scalefactor)
ax1.set_xlabel('Horizontal position $x$ (m)',fontsize=axistitlesize*scalefactor)
ax1.set_ylabel('Bed elevation $z_b$ (m)',fontsize=axistitlesize*scalefactor)
ax2.set_xlabel('Horizontal position $x$ (m)',fontsize=axistitlesize*scalefactor)
ax2.set_ylabel('Bed elevation $z_b$ (m)',fontsize=axistitlesize*scalefactor)
ax3.set_xlabel('Horizontal position $x$ (m)',fontsize=axistitlesize*scalefactor)
ax3.set_ylabel('Bed elevation $z_b$ (m)',fontsize=axistitlesize*scalefactor)

big_subplot.grid(visible=True, axis='y', color='lightgray', linewidth=0.1)
big_subplot.grid(visible=True, axis='x', color='lightgray', linewidth=0.1)

#Adds a legend to the plot
fig.legend(fontsize = legendsize * scalefactor, loc='upper center', title='                                                                      Bed elevation $z_b$ after running flood of hyperconcentration $c_{hyp}$ (%):',title_fontsize= legendtitlesize * scalefactor,frameon=False, ncol=4,bbox_to_anchor=(0.5, 1.0))#, bbox_to_anchor=(0.5, -0.03))

#Adjusts the plot layout
fig.tight_layout()
fig.subplots_adjust(top=0.92)

#Saves the figure
fig.savefig(save_file_path, dpi=1000)

#fig.show()






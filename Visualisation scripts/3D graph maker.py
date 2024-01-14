import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
from cycler import cycler
import sys
from mpl_toolkits.mplot3d import Axes3D
from colour import Color
plt.ioff()

labelpadsize = 0

plt.rcdefaults()
plt.rcParams["font.family"] = "Arial"
plt.rcParams['lines.linewidth'] = 0.1
plt.rcParams['axes.linewidth'] = 0.2
plt.rcParams['figure.dpi'] = 800
plt.rcParams['figure.figsize'] = 1.75*1.5,1.5*1.5
plt.rcParams['axes.prop_cycle'] = cycler(color=['blue','orange','brown','pink','gray','green','red','purple','#4B0082', '#2E8B57', '#7CFC00', '#C71585', '#808000', 
                                                '#228B22', '#008B8B', '#9ACD32', '#FF69B4', '#FF00FF', '#8B4513', '#708090', '#FFD700', '#006400', '#FF7F50', '#0000FF', 
                                                '#DC143C', '#1E90FF', '#FF0000', '#4682B4', '#556B2F', '#FFA07A', '#483D8B', '#BDB76B', '#FF4500', '#8A2BE2', '#CD5C5C', 
                                                '#8FBC8F', '#000080', '#FFFF00', '#8B0000', '#00FA9A', '#FF8C00'])

plt.rcParams['xtick.major.pad'] = -6
plt.rcParams['ytick.major.pad'] = -6
plt.rcParams['xtick.major.width'] = 0.2  # X-axis major tick width
plt.rcParams['ytick.major.width'] = 0.2
plt.rcParams['axes.titlepad'] = 0
plt.rcParams['axes.labelpad'] = -15

#Set to 1 for depth/bed_slope/free surface
#Set to 2 for single variable
#Set to 3 for erosion and deposition
display_case =  2

#Sets the individual variable for display case 2
#Set to 1 for Rouse number, 2 for concentration, 3 for bed elevation, 4 for entrainment
#5 for bed flowrate, 7 for deposition, 8 for change in bed elevation, 9 for median diameter
#10 for flowrate, 11 for water density, 12 for Courant number, 13 for velocity
# 14 for depth, 15 for flowrate
display_case_var = 3
plot_surface = 0
hide_depth = 'no'
depth_only = 'no'

starting_time = 0
finishing_time = 50000
line_rate = 10000
skip_rate = 50

extra_black_lines = 'no'
black_line_rate = 25

#Sets editable parameters for all charts
legendsize = 10
axistitlesize = 10
ticksize = 8
scalefactor = 0.25
titlesize = 15
axesspinecolour = 'gray'
allaxesspines = ['bottom','left','right','top']

#Defines the filepaths to acccess
filepath1 = r'../Coupled water-sediment model/output_depth.csv'
filepath2 = r'../Coupled water-sediment model/output_free_surface_elevation.csv'
filepath3 = r'../Coupled water-sediment model/output_bed_elevation.csv'
filepath4 = r'../Coupled water-sediment model/output_Rouse_number.csv'
filepath5 = r'../Coupled water-sediment model/output_concentration.csv'
filepath6 = r'../Coupled water-sediment model/output_entrainment.csv'
filepath7 = r'../Coupled water-sediment model/output_deposition.csv'
filepath8 = r'../Coupled water-sediment model/output_bed_flowrate.csv'
filepath9 = r'../Coupled water-sediment model/output_slope_diffusion.csv'
filepath10 = r'../Coupled water-sediment model/output_bed_elevation_change.csv'
filepath11 = r'../Coupled water-sediment model/output_median_diameter.csv'
filepath12 = r'../Coupled water-sediment model/output_flowrate.csv'
filepath13 = r'../Coupled water-sediment model/output_water_density.csv'
filepath14 = r'../Coupled water-sediment model/output_courant_number.csv'
filepath15 = r'../Coupled water-sediment model/output_velocity.csv'
filepath16 = r'../Coupled water-sediment model/output_flowrate.csv'


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

if display_case_1 == 1:

    filepath3 = r"C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 1, Sandbank\No_limiter_bed_elevation.csv"
    df_data = pd.read_csv(filepath3)
    n_timesteps = 50000
    n_steps = 100
    

#Converts the final time into the index of the row corresponding to that time
finishing_entry = int((n_timesteps/final_time) * finishing_time)
starting_entry = int((n_timesteps/final_time) * starting_time)



colour1 = Color("powderblue")
colour2 = Color("brown")
colour3 = Color("yellow")
colour4 = Color("green")
colour5 = Color("lightcyan")
colour6 = Color("bisque")
colour7 = Color("skyblue")

if starting_time == 0:

    colours = list(colour1.range_to(Color("darkviolet"),finishing_entry))
    colours2 = list(colour2.range_to(Color("red"),finishing_entry))
    colours3 = list(colour3.range_to(Color("brown"),finishing_entry))
    colours4 = list(colour4.range_to(Color("pink"),finishing_entry))
    colours5 = list(colour3.range_to(Color("darkgreen"),finishing_entry))
    colours6 = list(colour5.range_to(Color("darkturquoise"),finishing_entry))
    colours7 = list(colour6.range_to(Color("sienna"),finishing_entry))
    colours8 = list(colour7.range_to(Color("navy"),finishing_entry))

else:

    colours = list(colour1.range_to(Color("darkviolet"),finishing_entry - starting_entry))
    colours2 = list(colour2.range_to(Color("red"),finishing_entry - starting_entry))
    colours3 = list(colour3.range_to(Color("brown"),finishing_entry - starting_entry))
    colours4 = list(colour4.range_to(Color("pink"),finishing_entry - starting_entry))
    colours5 = list(colour3.range_to(Color("darkgreen"),finishing_entry - starting_entry))
    colours6 = list(colour5.range_to(Color("darkturquoise"),finishing_entry - starting_entry))
    colours7 = list(colour6.range_to(Color("sienna"),finishing_entry - starting_entry))
    colours8 = list(colour7.range_to(Color("navy"),finishing_entry - starting_entry))

matplotlib_colors = [c.rgb for c in colours] 
matplotlib_colors2 = [c.rgb for c in colours2] 
matplotlib_colors3 = [c.rgb for c in colours3] 
matplotlib_colors4 = [c.rgb for c in colours4] 
matplotlib_colors5 = [c.rgb for c in colours5] 
matplotlib_colors6 = [c.rgb for c in colours6] 
matplotlib_colors7 = [c.rgb for c in colours7] 
matplotlib_colors8 = [c.rgb for c in colours8]

 # Use z-values as colors
  # Choose a colormap (you can change it to any other colormap)

#Creates the figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') 

if parameter_case == 1:
    
    hide_surface = 'no'

if parameter_case == 2:
    
    hide_surface = 'no'
    
if parameter_case == 3:
    
    hide_surface = 'no'
    
if parameter_case == 4:
    
    hide_surface = 'yes'
    
if parameter_case == 5:
    
    hide_surface = 'yes'
    
if parameter_case == 6:
    
    hide_surface = 'yes'
    
if parameter_case == 7:
    
    hide_surface = 'yes'
    
if parameter_case == 8:
    
    hide_surface = 'yes'
    
if parameter_case == 9:
    
    hide_surface = 'yes'
    
if parameter_case == 10:
    
    hide_surface = 'yes'
    
if display_case == 1:

    #Reads the water depth output file to a dataframe
    df1 = pd.read_csv(filepath1,header=None,usecols=range(n_steps+5))
    #Reads the free surface elevation output file to a dataframe
    df2 = pd.read_csv(filepath2,header=None,usecols=range(n_steps+5))
    #Reads the bed elevation output file to a dataframe
    df3 = pd.read_csv(filepath3,header=None,usecols=range(n_steps+5))
    
    #For the water depth
    #The first row, corresponding to the x-positions
    y = df1.iloc[0, 1:]
    #The first column, corresponding to the time
    x = df1.iloc[1:, 0]
    #Each subsequent row, corresponding to the y-positions
    z = df1.iloc[1:, 1:].values  # Selecting rows excluding first and columns from 3rd onward
    
    #For the free surface elevation
    #The first row, corresponding to the x-positions
    y2 = df2.iloc[0, 1:]
    #The first column, corresponding to the time
    x2 = df2.iloc[1:, 0]
    #Each subsequent row, corresponding to the y-positions
    z2 = df2.iloc[1:, 1:].values  # Selecting rows excluding first and columns from 3rd onward
    
    #For the bed elevation
    #The first row, corresponding to the x-positions
    y3 = df3.iloc[0, 1:]
    #The first column, corresponding to the time
    x3 = df3.iloc[1:, 0]
    #Each subsequent row, corresponding to the y-positions
    z3 = df3.iloc[1:, 1:].values  # Selecting rows excluding first and columns from 3rd onward
    
    #Creates the 3D surface plot
    for i in range(starting_entry,finishing_entry):
        
        if starting_time == 0:
            color_pos = i
        else: 
            color_pos = int(i - starting_entry)
        try:
            if i % skip_rate == 0:
                if depth_only == 'no':
                    if i % line_rate == 0:
                        ax.plot([x.iloc[i]] * len(y),y3,z3[i],color = 'black',linewidth = 0.2)
                        ax.plot([x.iloc[i]] * len(y),y2,z2[i],color = 'black',linewidth = 0.2,zorder=2.2)
                        if hide_surface == 'no':
                            ax.plot([x.iloc[i]] * len(y),y,z[i],color = 'black',linewidth = 0.2)
                    else:
                        ax.plot([x.iloc[i]] * len(y),y3,z3[i],color=matplotlib_colors7[color_pos],linewidth = 0.1,zorder=0.9)
                        ax.plot([x.iloc[i]] * len(y),y2,z2[i],color=matplotlib_colors6[color_pos],linewidth = 0.1,zorder=2.1)
                        
                        if hide_surface == 'no':
                            ax.plot([x.iloc[i]] * len(y),y,z[i],color=matplotlib_colors8[i],linewidth = 0.1)
                            
                else:
                    if i % line_rate == 0:
                        ax.plot([x.iloc[i]] * len(y),y,z[i],color = 'black',linewidth = 0.2)
                    else:
                        ax.plot([x.iloc[i]] * len(y),y,z[i],color=matplotlib_colors7[i],linewidth = 0.1)
                        
            else:
                pass
            
        except ValueError:
            
            break
        
        except TypeError:
            
            break

    

    #Sets the axes limits depending on the case being run
    if parameter_case == 2 and display_case ==  1:  
        
        time = x.iloc[starting_entry:finishing_entry + 1].tolist()
        
        ax.set_zlim(bottom = 0,top=7.0)
        
        w_settling = 0.319106688674029
        
        conc_0 = 0.05
        
        time_constant = 5/(w_settling * (1 - conc_0/(1- 0.4)))

        for i in range(starting_entry,finishing_entry):
            
            
            t = x.iloc[i]
            
            operation_h = lambda x:  5 + (conc_0 * 5)/(1- 0.4 - conc_0) * (math.exp(-x/time_constant) -1)
            operation_z = lambda x:  1 - (conc_0 * 5)/(1- 0.4 - conc_0) * (math.exp(-x/time_constant) -1)
            
            h_analytical = [operation_h(item) for item in time]
            z_b_analytical = [operation_z(item) for item in time]
            
            #h_analytical = 5 - (conc_0 * 5)/(1- 0.4 - conc_0) * (math.exp(-t/time_constant -1))
            #z_b_analytical = 1 + (conc_0 * 5)/(1- 0.4 - conc_0) * (math.exp(-t/time_constant -1))
            
        #print(h_analytical,z_b_analytical)
        
        #print(total_length)

         
        for i in [250,500,750]:
            
            ax.plot(time,[i] * len(time),h_analytical,color = 'darkred',linewidth = 0.2,zorder=2,linestyle='dashed')
            ax.plot(time,[i] * len(time),z_b_analytical,color = 'darkred',linewidth = 0.2,zorder=0.95,linestyle='dashed')
            
            
            #sys.exit()
    
    if parameter_case == 3:
        
        time = x.iloc[starting_entry:finishing_entry + 1].tolist()

        
        operation_h = lambda x: 6 + 0.04/(1-0.4) * x
        operation_z = lambda x: 1 - 0.04/(1-0.4) * x
        
        h_analytical = [operation_h(item) for item in time]
        z_analytical = [operation_z(item) for item in time]
        
        
        for i in [250,500,750]:
            
            ax.plot(time,[i] * len(time),h_analytical,color = 'darkred',linewidth = 0.2,zorder=2,linestyle='dashed')
            ax.plot(time,[i] * len(time),z_analytical,color = 'darkred',linewidth = 0.2,zorder=0.95,linestyle='dashed')
            
    
    if parameter_case == 4:
        ax.set_zlim(bottom = 0,top=70)  
        
        if depth_only == 'yes':
            ax.set_zlim(bottom = 0,top=4)  
    
    ax.set_zlabel('Elevation (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
    
    if depth_only == 'yes':
        ax.set_zlabel('Flow depth (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
    
    chart_title = 'Temporal evolution of river'
    
if display_case == 2:
    
    if display_case_var == 1:
    
        filepath = filepath4
        ax.set_zlabel('Rouse number (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
    
    elif display_case_var == 2:
    
        filepath = filepath5
        ax.set_zlabel('Concentration (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 3:
    
        filepath = filepath3
        ax.set_zlabel('Bed elevation (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)

    elif display_case_var == 4:
    
        filepath = filepath6
        ax.set_zlabel('Entrainment (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 5:
    
        filepath = filepath8
        ax.set_zlabel('Bed flowrate (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 6:
    
        filepath = filepath8
        ax.set_zlabel('Slope diffusion (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
    
    elif display_case_var == 7:
    
        filepath = filepath7
        ax.set_zlabel('Deposition (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)

    elif display_case_var == 8:
    
        filepath = filepath10
        ax.set_zlabel('Change in bed elevation (m)', fontsize=axistitlesize*scalefactor,labelpad = labelpadsize)
        
    elif display_case_var == 9:
    
        filepath = filepath11
        ax.set_zlabel('Median diameter (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 10:
    
        filepath = filepath12
        ax.set_zlabel('Flowrate (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 11:
    
        filepath = filepath13
        ax.set_zlabel('Water density (kg m$-3$)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 12:
    
        filepath = filepath14
        ax.set_zlabel('Courant number', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 13:
    
        filepath = filepath15
        ax.set_zlabel('Velocity (ms$^-$$^1$)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)

    elif display_case_var == 14:
    
        filepath = filepath1
        ax.set_zlabel('Depth (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
        
    elif display_case_var == 15:
    
        filepath = filepath16
        ax.set_zlabel('Flowrate (m$^2$s$^-$$^1$)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
    

    #Reads the water depth output file to a dataframe
    df1 = pd.read_csv(filepath,header=None,usecols=range(n_steps+5))
    
    if display_case_var == 12:
        
        df1 = pd.read_csv(filepath,header=None,usecols=range(n_steps))

    
    #Creates the figure and 3D axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    #For the water depth
    #The first row, corresponding to the x-positions
    y = df1.iloc[0, 1:]
    #The first column, corresponding to the time
    x = df1.iloc[1:, 0]
    #Each subsequent row, corresponding to the y-positions
    z = df1.iloc[1:, 1:].values  # Selecting rows excluding first and columns from 3rd onward
    
    colormap = 'viridis'
    cmap = plt.get_cmap(colormap)

    x_grid, y_grid = np.meshgrid(x, y)
    
    if plot_surface == 1:
        surf = ax.plot_surface(x_grid.T, y_grid.T, z, cmap='viridis',alpha=0.5)    

    #Creates the 3D surface plot
    for i in range(starting_entry,finishing_entry):
        
        if starting_time == 0:
            color_pos = i
        else: 
            color_pos = int(i - starting_entry)
        
        try:
            if i % skip_rate == 0:
                
                    if i % line_rate == 0:
                        ax.plot([x.iloc[i]] * len(y),y,z[i],color = 'black',linewidth = 0.2, zorder = 20)
                    else:
                        if plot_surface != 1:
                            ax.plot([x.iloc[i]] * len(y),y,z[i],color=matplotlib_colors3[color_pos],linewidth = 0.1,alpha=0.4)
                        continue
                        
                       
                    
                    if extra_black_lines == 'yes':
                        #Mark the nth point in the y-z plane for each x-point in black
                        for j in range(2,len(y)):
                            
                            #print(j)
                            
                            if j % black_line_rate == 0:  # Change 'nth_point' to specify which point to highlight
                                
                                ax.scatter([x.iloc[i]], y[j], z[i], color='black',s=0.01,marker='x')  # Set the nth point in y-z plane to black
                        
                        
            else:
                pass
            
        except ValueError:
            
            break

    #Sets the axes limits depending on the case being run
    if parameter_case == 2 and display_case_var == 3:
        ax.set_zlabel('Volumetric sediment concentration', fontsize=axistitlesize*scalefactor)#,labelpad = 0.5)
        
    #if parameter_case == 5:
        #ax.set_zlim(bottom = 0,top=10) 
        
    if parameter_case == 1 and display_case_var == 6 :
            ax.set_zlim(bottom = 0,top=0.04)
            
    if parameter_case == 1 and display_case_var == 3:
        
        ax.set_zlim(bottom = 0,top=1.2) 
        ax.set_zlabel('Bed elevation $z_b$ (m)', fontsize=axistitlesize*scalefactor)
         

    chart_title = 'Temporal evolution of river'
    
    if plot_surface == 1:
        cbar = fig.colorbar(surf, ax=ax, orientation='horizontal', pad=-0.15, shrink=0.2)
        cbar.set_label('Change in bed elevation (m)')  # Add label to the color bar
        cbar.ax.xaxis.label.set_size(axistitlesize*scalefactor)
        cbar.ax.tick_params(labelsize = ticksize*scalefactor,length=2, pad=2, labeltop=False)
    

#The case for displaying the entrainment and deposition coefficients
if display_case == 3:

    #Reads the water depth output file to a dataframe
    df1 = pd.read_csv(filepath6,header=None,usecols=range(n_steps+5))
    #Reads the free surface elevation output file to a dataframe
    df2 = pd.read_csv(filepath7,header=None,usecols=range(n_steps+5))
    
    #For the water depth
    #The first row, corresponding to the x-positions
    y = df1.iloc[0, 1:]
    #The first column, corresponding to the time
    x = df1.iloc[1:, 0]
    #Each subsequent row, corresponding to the y-positions
    z = df1.iloc[1:, 1:].values  # Selecting rows excluding first and columns from 3rd onward
    
    #For the free surface elevation
    #The first row, corresponding to the x-positions
    y2 = df2.iloc[0, 1:]
    #The first column, corresponding to the time
    x2 = df2.iloc[1:, 0]
    #Each subsequent row, corresponding to the y-positions
    z2 = df2.iloc[1:, 1:].values  # Selecting rows excluding first and columns from 3rd onward
    
    #Creates the 3D surface plot
    for i in range(starting_entry,finishing_entry):
        
        try:
            if i % skip_rate == 0:
                
                    if i % line_rate != 0:

                        ax.plot([x.iloc[i]] * len(y),y2,z2[i],color = 'green',linewidth = 0.1)
                        ax.plot([x.iloc[i]] * len(y),y,z[i],color = 'red',linewidth = 0.1)
                        
                    else:
  
                        ax.plot([x.iloc[i]] * len(y),y2,z2[i],color = 'black',linewidth = 0.3)
                        ax.plot([x.iloc[i]] * len(y),y,z[i],color = 'black',linewidth = 0.3)   
                        
            else:
                pass
            
        except ValueError:
            
            break

    #Sets the axes limits depending on the case being run

    
    if parameter_case == 2:  
        ax.set_zlim(bottom = 0,top=7)
        
    if parameter_case == 4:
        ax.set_zlim(bottom = 0,top=70)  
        
    
    ax.set_zlabel('Elevation (m)', fontsize=axistitlesize*scalefactor)#,labelpad = labelpadsize)
    chart_title = 'Temporal evolution of river'

# Setting labels and title

ax.set_xlim(starting_time,right= finishing_time)  
ax.set_ylim(bottom = 0,top = total_length)

ax.xaxis._axinfo["grid"].update({"linewidth":0.1})
ax.yaxis._axinfo["grid"].update({"linewidth":0.1})
ax.zaxis._axinfo["grid"].update({"linewidth":0.1})

ax.tick_params(axis='both', labelsize=ticksize*scalefactor,width = 0.1)
#ax.xaxis.set_tick_params(width=0.2)
#ax.yaxis.set_tick_params(width=0.2)
ax.tick_params(axis='y', which='major', labelsize=ticksize*scalefactor,width = 0.1)
ax.tick_params(axis='x', which='major', labelsize=ticksize*scalefactor,width = 0.1)


ax.set_xlabel('Time (s)', fontsize=axistitlesize*scalefactor)#,labelpad = 0.5)
ax.set_ylabel('X-position (m)', fontsize=axistitlesize*scalefactor)#,labelpad = 0.5)

#ax.set_ylim(bottom =25000 ,top=35000)

if display_case_var == 8:

    filepath = filepath10
    ax.set_zlabel('Change in bed elevation (m)', fontsize=axistitlesize*scalefactor,labelpad = -15)
    
if display_case_var == 3:

    filepath = filepath10
    ax.set_zlabel('Bed elevation $z_b$ (m)', fontsize=axistitlesize*scalefactor,labelpad = -15)

#ax.set_title(chart_title,fontsize=titlesize*scalefactor,y=1)

ax.spines['bottom'].set_position(('axes', 0))  # Bottom spine at y=2
ax.spines['left'].set_position(('axes', 0))    # Left spine at x=1
ax.spines['right'].set_position(('axes', 0))
ax.spines['top'].set_position(('axes', 0))      # Left spine at x=1
ax.spines['top'].set_visible(True)           # Hiding the top spine
ax.spines['right'].set_visible(True)

#Sets the aspect ratio for the 3D graph
ax.set_box_aspect([1, 1, 0.26])

#Sets the graph view orientation
ax.view_init(elev=20, azim=-35) 

if parameter_case == 2:
    
    ax.set_box_aspect([1, 0.5, 0.5])
    ax.view_init(elev=15, azim=-60) 
    
if parameter_case == 3:
    
    ax.set_box_aspect([1, 0.5, 0.5])
    ax.view_init(elev=15, azim=-60) 
    
#fig.tight_layout()
# Show plot
fig.show()

save_file_path = r'C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Diagrams\Graphical outputs\Case 1, Sandbank\Bed_elevation_failure.png'  # Replace this with your desired file path
png_name = input("Enter the file name of the png to save: ")
save_file_path = r"../Results/Graphical results/{}.png".format(png_name)
plt.savefig(save_file_path, dpi=1000)

print('Figure saved!')
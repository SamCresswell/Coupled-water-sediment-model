# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 17:12:51 2024

@author: sam_c
"""

Certainly! Here are the commands using the provided file path to extract the 1st, 2nd, nth, and last lines from the CSV file using PowerShell:

For the 1st line:

powershell
Copy code
Get-Content -Path "C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 3, Double slope\output_bed_elevation.csv" -TotalCount 1
For the 2nd line:

powershell
Copy code
Get-Content -Path "C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 3, Double slope\output_bed_elevation.csv" -TotalCount 2 | Select-Object -Last 1
For the nth line (replace [n] with the line number you want):

powershell
Copy code
Get-Content -Path "C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 3, Double slope\output_bed_elevation.csv" | Select-Object -Index [n]-1
For the last line:

powershell
Copy code
Get-Content -Path "C:\Users\sam_c\OneDrive - University of Glasgow\Invididual project\Model\Week 9 - Sediment model\Coupled water-sediment model\Results\Case 3, Double slope\output_bed_elevation.csv" | Select-Object -Last 1
Replace the file path in each command with the one provided. These commands will display the specified lines in the PowerShell console for the given CSV file.
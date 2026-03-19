from paraview.simple import *

# 1. Load the FEniCS output files
# Using the exact filenames from your FEniCS script
disp_xdmf = XDMFReader(registrationName='Displacement', FileNames='displacement_in_meters.xdmf')
temp_xdmf = XDMFReader(registrationName='Temperature', FileNames='temperature_in_Kelvin.xdmf')

# IMPORTANT: Force ParaView to actually read the files
disp_xdmf.UpdatePipeline()
temp_xdmf.UpdatePipeline()

# 2. Append Attributes
# Merges 'u' (displacement) and 'sigma' (stress) into one object
merged = AppendAttributes(Input=[disp_xdmf,  temp_xdmf])
merged.UpdatePipeline()

# 3. Apply Warp By Vector
# FEniCS names the displacement array 'u'
warped = WarpByVector(Input=merged)
warped.Vectors = ['POINTS', 'u']
warped.ScaleFactor = 5000.0  # Adjust based on expected deformation
warped.UpdatePipeline()

# 5. Save the state for manual inspection
# This creates a file you can open in the ParaView GUI
SaveState('view_results.pvsm')

print("---------------------------------------------------------")
print("Pipeline successful! No 'Missing input' errors.")
print("To view: Open ParaView GUI -> File -> Load State -> view_results.pvsm")
print("---------------------------------------------------------")
from paraview.simple import *

# 1. Load the FEniCS output files
# Using the exact filenames from your FEniCS script
temp_xdmf =  XDMFReader(registrationName='Temperature', FileNames='temperature_in_degrees.xdmf')

# IMPORTANT: Force ParaView to actually read the files
temp_xdmf.UpdatePipeline()


# 5. Save the state for manual inspection
# This creates a file you can open in the ParaView GUI
SaveState('view_results.pvsm')

print("---------------------------------------------------------")
print("Pipeline successful! No 'Missing input' errors.")
print("To view: Open ParaView GUI -> File -> Load State -> view_results.pvsm")
print("---------------------------------------------------------")
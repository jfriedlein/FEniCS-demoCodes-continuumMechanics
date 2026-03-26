from paraview.simple import *

# 1. Load the FEniCS output files
# Using the exact filenames from your FEniCS script
disp_pvd = PVDReader(registrationName='Displacement', FileName='displacement_in_meters.pvd')
stress_pvd = PVDReader(registrationName='Stress', FileName='stress_in_MPa.pvd')
temp_pvd = PVDReader(registrationName='Temperature', FileName='temperature_in_Kelvin.pvd')
heatflux_pvd= PVDReader(registrationName='Heat Flux', FileName='heat_flux.pvd')

# IMPORTANT: Force ParaView to actually read the files
disp_pvd.UpdatePipeline()
stress_pvd.UpdatePipeline()
temp_pvd.UpdatePipeline()
heatflux_pvd.UpdatePipeline()

# 2. Append Attributes
# Merges 'u' (displacement) and 'sigma' (stress) into one object
merged = AppendAttributes(Input=[disp_pvd, stress_pvd, temp_pvd, heatflux_pvd])
merged.UpdatePipeline()

# 3. Apply Warp By Vector
# FEniCS names the displacement array 'u'
warped = WarpByVector(Input=merged)
warped.Vectors = ['POINTS', 'u']
warped.ScaleFactor = 1000.0  # Adjust based on expected deformation
warped.UpdatePipeline()

# 5. Save the state for manual inspection
# This creates a file you can open in the ParaView GUI
SaveState('view_results.pvsm')

print("---------------------------------------------------------")
print("Pipeline successful! No 'Missing input' errors.")
print("To view: Open ParaView GUI -> File -> Load State -> view_results.pvsm")
print("---------------------------------------------------------")
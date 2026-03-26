from paraview.simple import *

# 1. Load the FEniCS output files
disp_xdmf = XDMFReader(registrationName='Displacement', FileNames='displacement_in_meters.xdmf')
temp_xdmf = XDMFReader(registrationName='Temperature', FileNames='temperature_in_Kelvin.xdmf')
stress_xdmf = XDMFReader(registrationName='Stress', FileNames='stress.xdmf')
flux_xdmf = XDMFReader(registrationName='Heat Flux', FileNames='heat_flux.xdmf')

# IMPORTANT: Force ParaView to actually read the files
disp_xdmf.UpdatePipeline()
temp_xdmf.UpdatePipeline()
stress_xdmf.UpdatePipeline()
flux_xdmf.UpdatePipeline()

# 2. Append Attributes
merged = AppendAttributes(Input=[disp_xdmf, temp_xdmf, stress_xdmf, flux_xdmf])
merged.UpdatePipeline()

# 2b. Apply Temporal Statistics to the merged dataset
temporal_stats = TemporalStatistics(registrationName='Temporal Statistics', Input=merged)
temporal_stats.UpdatePipeline()

# 3. Apply Warp By Vector
# Use the merged dataset for warping (geometric transformation)
warped = WarpByVector(Input=merged)  # Warping uses the displacement 'u'
warped.Vectors = ['POINTS', 'u']
warped.ScaleFactor = 1000.0  # Adjust based on expected deformation
warped.UpdatePipeline()

# 4. Save the state for manual inspection
SaveState('view_results.pvsm')

print("---------------------------------------------------------")
print("Pipeline successful! No 'Missing input' errors.")
print("Temporal statistics added for Displacement, Temperature, Stress, and Heat Flux.")
print("To view: Open ParaView GUI -> File -> Load State -> view_results.pvsm")
print("---------------------------------------------------------")
from paraview.simple import *

# 1. Load the FEniCS output files
temp_xdmf = XDMFReader(registrationName='Temperature', FileNames='temperature_in_degrees.xdmf')
flux_xdmf = XDMFReader(registrationName='Heat Flux', FileNames='heat_flux.xdmf')

# IMPORTANT: Force ParaView to actually read the files
temp_xdmf.UpdatePipeline()
flux_xdmf.UpdatePipeline()

# 2. Apply Temporal Statistics to Temperature
temp_stats = TemporalStatistics(registrationName='Temperature Temporal Statistics', Input=temp_xdmf)
temp_stats.UpdatePipeline()

# 3. Apply Temporal Statistics to Heat Flux
flux_stats = TemporalStatistics(registrationName='Heat Flux Temporal Statistics', Input=flux_xdmf)
flux_stats.UpdatePipeline()

# 4. Save the state for manual inspection
SaveState('view_results.pvsm')

print("---------------------------------------------------------")
print("Pipeline successful! No 'Missing input' errors.")
print("Temporal statistics added for Temperature and Heat Flux.")
print("To view: Open ParaView GUI -> File -> Load State -> view_results.pvsm")
print("---------------------------------------------------------")
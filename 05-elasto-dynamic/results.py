from paraview.simple import *
import os
import sys

# Names of the XDMF files
disp_filename = 'displacement_in_meters.xdmf'
stress_filename = 'stress_in_MPa.xdmf'

# First, check current directory
current_folder = '.'
disp_file = os.path.join(current_folder, disp_filename)
stress_file = os.path.join(current_folder, stress_filename)

# Fallback folder
fallback_folder = 'ev0'

if not (os.path.isfile(disp_file) and os.path.isfile(stress_file)):
    print(f"Files not found in current directory. Trying fallback folder '{fallback_folder}'...")
    disp_file = os.path.join(fallback_folder, disp_filename)
    stress_file = os.path.join(fallback_folder, stress_filename)
    if not (os.path.isfile(disp_file) and os.path.isfile(stress_file)):
        print(f"Error: Required XDMF files not found in current directory or '{fallback_folder}'.")
        sys.exit(1)

print(f"Using displacement file: {disp_file}")
print(f"Using stress file: {stress_file}")

# 1. Load the FEniCS output files
disp_xdmf = XDMFReader(registrationName='Displacement', FileNames=[disp_file])
stress_xdmf = XDMFReader(registrationName='Stress', FileNames=[stress_file])

# IMPORTANT: Force ParaView to actually read the files
disp_xdmf.UpdatePipeline()
stress_xdmf.UpdatePipeline()

# 2. Append Attributes
merged = AppendAttributes(Input=[disp_xdmf, stress_xdmf])
merged.UpdatePipeline()

# 2b. Apply Temporal Statistics
temporal_stats = TemporalStatistics(registrationName='Temporal Statistics', Input=merged)
temporal_stats.UpdatePipeline()

# 3. Apply Warp By Vector
warped = WarpByVector(Input=merged)  # Use 'merged', not temporal_stats, to warp the geometry
warped.Vectors = ['POINTS', 'u']
warped.ScaleFactor = 1500.0
warped.UpdatePipeline()

# 4. Save the state for manual inspection
SaveState('view_results.pvsm')

print("---------------------------------------------------------")
print("Pipeline successful! No 'Missing input' errors.")
print("Temporal statistics added for Displacement and Stress.")
print("To view: Open ParaView GUI -> File -> Load State -> view_results.pvsm")
print("---------------------------------------------------------")
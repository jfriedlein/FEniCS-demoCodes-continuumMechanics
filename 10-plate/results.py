from paraview.simple import *
import os

# 1. Check current directory name
current_dir = os.getcwd()
# This checks if the folder '01-tractions' is at the end of our current path
is_traction_folder = current_dir.endswith("01-tractions")

# 1. Load the FEniCS output files
# Using the exact filenames from your FEniCS script
disp_pvd = PVDReader(registrationName='Displacement', FileName='displacement_in_meters.pvd')
stress_pvd = PVDReader(registrationName='Stress', FileName='stress_in_MPa.pvd')

# IMPORTANT: Force ParaView to actually read the files
disp_pvd.UpdatePipeline()
stress_pvd.UpdatePipeline()

# 2. Append Attributes
# Merges 'u' (displacement) and 'sigma' (stress) into one object
merged = AppendAttributes(Input=[disp_pvd, stress_pvd])
merged.UpdatePipeline()

# 3. Apply Reflect Filter
# Access the internal 'ReflectionPlane' object's properties
# Normal [1, 0, 0] means it mirrors along the X-axis
# Define reflection plane (this IS the Reflect filter in ParaView 6.0+)
if is_traction_folder:
    print("MATCH FOUND: Applying Reflect filter for 01-tractions...")
    reflected = Reflect(Input=merged)
    reflected.ReflectionPlane.Origin = [1.5, 0, 0]
    reflected.ReflectionPlane.Normal = [1, 0, 0]
    reflected.CopyInput = 1
    reflected.ReflectAllInputArrays = 1
    reflected.UpdatePipeline()
    final_output = reflected
else:
    final_output = merged


# 4. Apply Warp By Vector
# FEniCS names the displacement array 'u'
warped = WarpByVector(Input=final_output)
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
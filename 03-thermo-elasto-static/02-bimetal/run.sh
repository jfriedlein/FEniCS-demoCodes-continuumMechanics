#!/bin/bash
# Setup and run the LKM demo in a local Conda environment (no Docker)

cd ..
cd ..

# 1. Create environment (skip if already exists)
cd environment/local
if ! conda info --envs | grep -q "fenics-lkm-env"; then
    conda env create -f environment.yml
else
    echo "Environment already exists, skipping creation."
fi

cd ..
cd ..

# 2. Run main.py
cd 03-thermo-elasto-static/02-bimetal

echo "Running the LKM demo..."
conda run -n fenics-lkm-env python main.py

if [ $? -eq 0 ]; then
    echo "Done! Check displacement.pvd, temperature.pvd and stress.pvd files."
else
    echo "ERROR: The LKM demo failed to run. Check the error messages above."
fi

# 3. Run results.py to generate ParaView state file
echo "Generating ParaView state file..."
conda run -n fenics-lkm-env pvbatch results.py
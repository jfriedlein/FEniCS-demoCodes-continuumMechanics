#!/bin/bash
# Setup and run the LKM demo

# 1. Create environment (skip if already exists)
if ! conda info --envs | grep -q "fenics-lkm-env"; then
    conda env create -f environment.yml
else
    echo "Environment already exists, skipping creation."
fi

# 2. Run main.py
echo "Running the LKM demo..."
conda run -n fenics-lkm-env python main.py

if [ $? -eq 0 ]; then
    echo "Done! Check displacement.pvd and stress.pvd files."
else
    echo "ERROR: The LKM demo failed to run. Check the error messages above."
fi
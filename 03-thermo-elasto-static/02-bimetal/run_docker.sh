#!/bin/bash
#setup and run the LKM demo in a Docker container

# Load Docker image
cd ..
cd ..

IMAGE_NAME="fenics-tool-offline:v1"

echo "[STEP 1] Loading Docker image if needed..."
docker image inspect $IMAGE_NAME >/dev/null 2>&1 || docker load -i environment/docker/fenics_tool.tar
echo "   Done"
echo ""

# Run FEniCS simulation
echo "[STEP 2] Running FEniCS simulation (displacement, stress & temperature)..."
cd 03-thermo-elasto-static/02-bimetal
docker run --rm -v $(pwd):/app fenics-tool-offline:v1 conda run -n fenics-lkm-env python /app/main.py
if [ $? -eq 0 ]; then echo "   Done"; else exit 1; fi
echo ""

# Process and visualize with ParaView
echo "[STEP 3] Processing results with ParaView..."
docker run --rm -v "$(pwd)":/app fenics-tool-offline:v1 conda run -n fenics-lkm-env pvbatch /app/results.py
if [ $? -eq 0 ]; then
    echo "   Workflow completed!"
    echo ""
    echo " Output: view_results.pvsm (open in ParaView)"
    echo ""
else
    exit 1
fi
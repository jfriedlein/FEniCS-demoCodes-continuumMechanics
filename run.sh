#!/bin/bash

# 1. Store the root directory
BASE_DIR=$(pwd)

echo "Current Root: $BASE_DIR"
echo "Enter the relative path to the folder (e.g., 01-elasto-static/01-tractions):"
read script_path

# 2. Define the absolute path to the target folder
TARGET_DIR="$BASE_DIR/$script_path"
# NEW: Define PARENT_DIR so we can find results.py outside the target
PARENT_DIR=$(dirname "$TARGET_DIR")

# Check if main.py exists in that folder
if [[ ! -f "$TARGET_DIR/main.py" ]]; then
    echo "ERROR: 'main.py' not found in $TARGET_DIR!"
    exit 1
fi

echo "Do you want to execute in Docker? [y/n]:"
read input

if [[ "$input" == "y" || "$input" == "Y" ]]; then
    
    IMAGE_NAME="fenics-tool-offline:v1"

    echo "[STEP 1] Loading Docker image if needed..."
    docker image inspect $IMAGE_NAME >/dev/null 2>&1 || docker load -i "$BASE_DIR/environment/docker/fenics_tool.tar"

    # [CRITICAL CHANGE]
    # We mount $TARGET_DIR directly to /app. 
    # This means Python's "current directory" inside Docker is your subfolder.
    
    echo "[STEP 2] Running Simulation (Files will save to $script_path)..."
    docker run --rm -v "$TARGET_DIR":/app $IMAGE_NAME conda run -n fenics-lkm-env python /app/main.py
    
    if [ $? -eq 0 ]; then
        echo "[STEP 3] Running ParaView results..."
        docker run --rm -v "$TARGET_DIR":/app -v "$PARENT_DIR":/scripts -w /app  $IMAGE_NAME conda run -n fenics-lkm-env pvbatch /scripts/results.py
        echo "Done! Files are saved in: $script_path"
    else 
        echo "Simulation failed."
        exit 1 
    fi

elif [[ "$input" == "n" || "$input" == "N" ]]; then
    # LOCAL EXECUTION
    if ! conda info --envs | grep -q "fenics-lkm-env"; then
        conda env create -f "$BASE_DIR/environment/local/environment.yml"
    fi

    # Move into the directory so Python saves files locally
    cd "$TARGET_DIR" || exit
    
    echo "Running locally in $(pwd)..."
    conda run -n fenics-lkm-env python main.py
    
    
    conda run -n fenics-lkm-env pvbatch $BASE_DIR/results.py
    
    
    echo "Done! Files saved in $(pwd)"
    
    # Return to root
    cd "$BASE_DIR"
fi
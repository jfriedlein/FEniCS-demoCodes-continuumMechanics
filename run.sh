#!/bin/bash

# 1. Store the root directory
BASE_DIR=$(pwd)

echo "Current Root: $BASE_DIR"
echo "Enter the relative path to the folder (e.g., 01-elasto-static/01-tractions):"
read -e script_path

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

    echo "[STEP 1/3] Loading Docker image if needed..."
    docker image inspect $IMAGE_NAME >/dev/null 2>&1 || docker load -i "$BASE_DIR/environment/docker/fenics_tool.tar"

    # [CRITICAL CHANGE]
    # We mount $TARGET_DIR directly to /app. 
    # This means Python's "current directory" inside Docker is your subfolder.
    
    echo "[STEP 2/3] Running Simulation (Files will save to $script_path)..."
    docker run --rm -v "$TARGET_DIR":/app $IMAGE_NAME conda run -n fenics-lkm-env python /app/main.py
    
    if [ $? -ne 0 ]; then
        echo "Simulation failed."
        exit 1
    fi

    #  choose results.py from subfolder first, then root
    if [[ -f "$TARGET_DIR/results.py" ]]; then
        RESULTS_LOCAL="$TARGET_DIR/results.py"
        RESULTS_CONTAINER="/app/results.py"
        MOUNT_RESULTS="-v $TARGET_DIR:/app"
    elif [[ -f "$PARENT_DIR/results.py" ]]; then
        RESULTS_LOCAL="$PARENT_DIR/results.py"
        RESULTS_CONTAINER="/scripts/results.py"
        MOUNT_RESULTS="-v $TARGET_DIR:/app -v $PARENT_DIR:/scripts"
    else
        echo "ERROR: results.py not found in $TARGET_DIR or parent."
        exit 1
    fi

    echo "[STEP 3/3] Running ParaView results..."
    docker run --rm $MOUNT_RESULTS -w /app $IMAGE_NAME conda run -n fenics-lkm-env pvbatch $RESULTS_CONTAINER
    echo "Done! Files are saved in: $script_path"

elif [[ "$input" == "n" || "$input" == "N" ]]; then
    # LOCAL EXECUTION
    echo "[STEP 1/3] creating local environment if not exists..."
    if ! conda info --envs | grep -q "fenics-lkm-env"; then
        conda env create -f "$BASE_DIR/environment/local/environment.yml"
    fi

    # Move into the directory so Python saves files locally
    cd "$TARGET_DIR" || exit
    
    echo "[STEP 2/3] Running locally in $(pwd)..."
    conda run -n fenics-lkm-env python main.py

    # Step 3: choose results.py from subfolder first, then root
    if [[ -f "$TARGET_DIR/results.py" ]]; then
        RESULTS_LOCAL="$TARGET_DIR/results.py"
    elif [[ -f "$PARENT_DIR/results.py" ]]; then
        RESULTS_LOCAL="$PARENT_DIR/results.py"
    else
        echo "ERROR: results.py not found in $TARGET_DIR or parent."
        cd "$BASE_DIR"
        exit 1
    fi

    echo "[STEP 3/3] Running ParaView results..."
    conda run -n fenics-lkm-env pvbatch "$RESULTS_LOCAL"

    echo "Done! Files saved in $(pwd)"

    # Return to root
    cd "$BASE_DIR"
fi
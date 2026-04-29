#!/bin/bash

BASE_DIR=$(pwd)

SCRIPT_PATH=$1
MODE=$2   # "docker" or "local"

if [[ -z "$SCRIPT_PATH" || -z "$MODE" ]]; then
    echo "Direct usage: ./run.sh <path> <docker|local>"
    echo "Starting command line usage ..."
    echo "Current Root: $BASE_DIR"
    echo "Enter the relative path to the folder (e.g., 01-elasto-static/01-tractions):"
    read -e SCRIPT_PATH
    echo "Do you want to execute in Docker? [y/n]:"
    read input
    if [[ "$input" == "y" ]]; then
        MODE="docker"
    else
        MODE="local"
    fi
fi

TARGET_DIR="$BASE_DIR/$SCRIPT_PATH"
PARENT_DIR=$(dirname "$TARGET_DIR")

if [[ ! -f "$TARGET_DIR/main.py" ]]; then
    echo "ERROR: main.py not found in $TARGET_DIR!"
    exit 1
fi

if [[ "$MODE" == "docker" ]]; then

    IMAGE_NAME="fenics-tool-offline:v1"

    echo "[STEP 1/3] Loading Docker image if needed..."
    docker image inspect $IMAGE_NAME >/dev/null 2>&1 || docker load -i "$BASE_DIR/environment/docker/fenics_tool.tar"

    echo "[STEP 2/3] Running Simulation..."
    docker run --rm -v "$TARGET_DIR":/app $IMAGE_NAME \
        conda run -n fenics-lkm-env python /app/main.py

    if [ $? -ne 0 ]; then
        echo "Simulation failed."
        exit 1
    fi

    if [[ -f "$TARGET_DIR/results.py" ]]; then
        RESULTS_CONTAINER="/app/results.py"
        MOUNT_RESULTS="-v $TARGET_DIR:/app"
    elif [[ -f "$PARENT_DIR/results.py" ]]; then
        RESULTS_CONTAINER="/scripts/results.py"
        MOUNT_RESULTS="-v $TARGET_DIR:/app -v $PARENT_DIR:/scripts"
    else
        echo "ERROR: results.py not found."
        exit 1
    fi

    echo "[STEP 3/3] Running ParaView..."
    docker run --rm $MOUNT_RESULTS -w /app $IMAGE_NAME \
        conda run -n fenics-lkm-env pvbatch $RESULTS_CONTAINER

    echo "Done! Files saved in: $SCRIPT_PATH"

else
    echo "[STEP 1/3] Checking environment..."
    if ! conda info --envs | grep -q "fenics-lkm-env"; then
        conda env create -f "$BASE_DIR/environment/local/environment.yml"
    fi

    cd "$TARGET_DIR" || exit

    echo "[STEP 2/3] Running locally..."
    conda run -n fenics-lkm-env python main.py

    if [[ -f "$TARGET_DIR/results.py" ]]; then
        RESULTS_LOCAL="$TARGET_DIR/results.py"
    elif [[ -f "$PARENT_DIR/results.py" ]]; then
        RESULTS_LOCAL="$PARENT_DIR/results.py"
    else
        echo "ERROR: results.py not found."
        cd "$BASE_DIR"
        exit 1
    fi

    echo "[STEP 3/3] Running ParaView..."
    conda run -n fenics-lkm-env pvbatch "$RESULTS_LOCAL"

    echo "Done! Files saved in $(pwd)"
    cd "$BASE_DIR"
fi

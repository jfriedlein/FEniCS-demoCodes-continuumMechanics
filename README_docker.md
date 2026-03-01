
## Use Docker for an isolated, reproducible environment without installing FEniCS locally.
# Linear Continuum Mechanics (LKM) Demo – FEniCS

This project demonstrates 3D linear elasticity using Python and FEniCS.

## Project Structure

├── main.py                         # FEniCS code 
├── results.py                      # ParaView visualization pipeline script
├── environment.yml                 # Conda environment definition
├── run_docker.sh                   # Docker execution script (Linux / WSL)
├── fenics_tool.tar                 # Docker container
├── displacement_in_meters.pvd      # Output: displacement field (generated)
├── stress_in_MPa.pvd               # Output: stress field (generated)
├── view_results.pvsm               # Output: ParaView state file (generated)
└── README_docker.md                # This documentation

### Prerequisites

- **if Windows**: Install WSL (Ubuntu)

## Docker Workflow

Run the complete workflow (simulation + visualization) inside Docker containers:

```bash
docker load -i fenics_tool.tar
```

### Run FEniCS Simulation
```bash
docker run --rm -v $(pwd):/app fenics-tool-offline:v1
```
Computes displacement and stress fields for the 3D linear elasticity problem.

### Process Results with ParaView
```bash
docker run --rm -v "$(pwd)":/app fenics-lkm-env conda run -nfenics-lkm-env pvbatch results.py
```
Creates reflected geometry, applies deformation, and generates `view_results.pvsm`.


## Quick Start

Use the provided script to run all steps automatically:
```bash
bash run_docker.sh
```

## Viewing Results

1. Open ParaView on your local machine
2. Go to: **File → Load State**
3. Select: **view_results.pvsm**
4. Examine displacement and stress visualizations




## Use Docker for an isolated, reproducible environment without installing FEniCS locally.
# Linear Continuum Mechanics (LKM) Demo – FEniCS

This project demonstrates 3D linear elasticity using Python and FEniCS.

## Project Structure

FENICS-LINEAR-CONTINUUM/
├── environment/
│   ├── docker/                 # REPRODUCIBILITY: Container configuration for running simulations in Docker
│   │   ├── environment.yml     # List of required dependencies (FEniCS, ParaView, Python packages)
│   │   └── fenics.tool.tar     # Prebuilt Docker image or instructions to build/load the container
│   └── local/                  # LOCAL SETUP: Configuration for running the project on a local machine
│       └── environment.yml     # Conda/Mamba environment specification for native installation
│
├── 01-elasto-static/           # PHYSICS MODULE: Linear elasticity simulations (static problems)
│   ├── 01-tractions/           # SCENARIO: Solid subjected to boundary traction forces
│   │   ├── main.py             # EXECUTION: Defines the FEM problem, boundary conditions, and solver routine
│   │   ├── results.py          # POST-PROCESSING: Extracts data and generates plots or derived quantities
│   │   ├── view_results.pvsm   # VISUALIZATION: Saved ParaView state for quick visualization of results
│   │   └── outputs/            # AUTO-GENERATED: Simulation output files
│   │       ├── paraview/       # Mesh and solution files (.pvd, .xdmf, .vtu) for visualization, generated after simulation.
│   │       └── plots/          # Generated figures (PNG/PDF) such as stress fields or convergence plots
│   │
│   └── 02-gravity/             # SCENARIO: Elastic body subjected to body forces (e.g., gravity/self-weight)
│       └── ...                 # Follows the same structure as other simulation scenarios
│
├── slides/                     # THEORY: Lecture or presentation material related to the simulations
│   └── lkm-demo.pdf            # Example slide deck describing LKM problems
│
├── README.md                   # DOCUMENTATION: Project overview, setup instructions, and usage guide
│
└── run.sh                      # EXECUTION SCRIPT: Helper script to run selected simulation modules

### Prerequisites
- **Install Docker**: (https://docs.docker.com/desktop/setup/install/linux/ubuntu/#install-docker-desktop)
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


## Local Setup Instructions

### Prerequisites

- **Linux / macOS**: Install Miniconda or Anaconda
- **Windows**: Install WSL (Ubuntu), then install Miniconda inside WSL

### Installation Steps

1. Create and activate the Conda environment:
```bash
conda env create -f environment.yml
conda activate fenics-lkm-env
```

2. Run the FEniCS simulation:
```bash
python main.py
```

3. Process results with ParaView:
```bash
conda install -c conda-forge paraview
pvbatch results.py
```


## Quick Start

Use the provided script to run all steps automatically:

```bash
bash run.sh
```

## Viewing Results

1. Open ParaView on your local machine
2. Go to: **File → Load State**
3. Select: **view_results.pvsm**
4. Examine displacement and stress visualizations



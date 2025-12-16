# Linear Continuum Mechanics (LKM) Demo – FEniCS

This project demonstrates a simple 2D linear elasticity example using Python and FEniCS.

## Project Structure

## Project Structure

```text
.
├── main.py              # FEniCS code for 2D linear elasticity
├── environment.yml      # Conda environment definition (minimal FEniCS setup)
├── run.sh               # One-command setup and execution script (Linux / WSL)
├── displacement.pvd     # Output: displacement field (generated)
├── stress.pvd           # Output: stress field (generated)
└── README.md            # Project documentation

```

## Setup Instructions

1. Install Conda
- Linux / macOS: Install Miniconda or Anaconda
- Windows: Install WSL (Ubuntu), then install Miniconda *inside WSL*

2. Create and activate the environment

Run the following commands on Linux, macOS, or WSL:

```bash
conda env create -f environment.yml
conda activate fenics-lkm-env
python main.py

```

3. Run the simulation
## Optional: One-Command Setup and Run
### Linux / WSL

```bash
chmod +x run.sh
./run.sh
```
3. View Results

# Results can be viewed using paraview

# Optional: install via conda if needed
conda install -n fenics-lkm-env -c conda-forge paraview  




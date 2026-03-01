# Linear Continuum Mechanics (LKM) Demo – FEniCS

This project demonstrates 3D linear elasticity using Python and FEniCS.

## Project Structure

```text
.
├── main.py                         # FEniCS code for 3D linear elasticity simulation
├── results.py                      # ParaView visualization pipeline script
├── environment.yml                 # Conda environment definition
├── run.sh                          # Local execution script (Linux / WSL)
├── displacement_in_meters.pvd      # Output: displacement field (generated)
├── stress_in_MPa.pvd               # Output: stress field (generated)
└── README.md                       # This documentation
```

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

## Quick Start (One Command)

### Linux / WSL
```bash
chmod +x run.sh
./run.sh
```

## Viewing Results

1. Open ParaView on your local machine
2. Go to: **File → Load State**
3. Select: **view_results.pvsm**
4. Examine displacement and stress visualizations






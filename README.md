# Linear Continuum Mechanics (LKM) Demo – FEniCS

This project provides a collection of 3D linear elasticity simulations using **FEniCS**, with support for both **local execution** and **Docker-based reproducible environments**.

It includes a **GUI-based launcher** and a **script-based workflow** to simplify running simulations and visualizing results.

---

# 📁 Project Structure

```
FENICS-LINEAR-CONTINUUM/
├── environment/
│   ├── docker/                 # Docker environment for reproducible runs
│   │   ├── environment.yml
│   │   └── fenics_tool.tar
│   └── local/                  # Local Conda environment
│       └── environment.yml
│
├── 01-elasto-static/
│   ├── 01-tractions/
│   │   ├── main.py
│   │   ├── results.py
│   │   ├── view_results.pvsm
│   │   └── outputs/
│   │
│   └── 02-gravity/
│       └── ...
│
├── 04-lshape/                  # Example WITHOUT subfolder (direct execution)
│   ├── main.py
│   ├── results.py
│   ├── *.vtu / *.pvd / *.png
│
├── gui.py                      # Graphical interface to browse and run simulations
├── run.sh                      # Main execution script (used by GUI and CLI)
│
└── README.md
```

---

# 🚀 How to Run

## Option 1 — GUI (Recommended)

Launch the graphical interface:

```bash
python exampleSelector.py
```

### Features:

* Browse topics and examples visually
* Supports both:

  * `topic/example/`
  * `topic/` (direct execution folders)
* Run simulations (Docker or local)
* Open results folder
* Launch ParaView (Windows automatically if using WSL)

---

## Option 2 — Script (CLI)

Run using:

```bash
bash run.sh <path> <docker|local>
```

### Examples:

```bash
bash run.sh 01-elasto-static/01-tractions docker
bash run.sh 04-lshape local
```

---

# Docker Workflow (Recommended for reproducibility)

### Load Docker image (first time only)

```bash
docker load -i environment/docker/fenics_tool.tar
```

### Run simulation + visualization

Handled automatically via:

```bash
bash run.sh <path> docker
```

---

# 💻 Local Execution

### Setup environment

```bash
conda env create -f environment/local/environment.yml
conda activate fenics-lkm-env
```

### Run manually by

```bash
cd <example-folder>
python main.py
pvbatch results.py
```

---

# 📊 Results & Visualization

After execution:

* Results are saved in:

  * `outputs/`, `results/`, or directly in the example folder
* Visualization file:

  * `view_results.pvsm`

---

## 🔍 Open in ParaView

### GUI:

Click **“Open ParaView”**

### Manual:

1. Open ParaView
2. File → Load State
3. Select `view_results.pvsm`

---

# 🪟 Windows + WSL Support

* GUI automatically:

  * Converts Linux paths → Windows paths
  * Opens folders in Windows Explorer
  * Launches ParaView on Windows

---


# 📦 Requirements

## Required

* Docker (recommended)
* OR Conda / Mamba

## Windows users

* Install **WSL (Ubuntu)**
* Install Docker Desktop

---

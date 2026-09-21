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

# How to Setup and Install
* Download the repository as .zip
* Unzip the downloaded file "FEniCS-demoCodes-continuumMechanics-master.zip"
* Download and install Docker following the instructions on https://docs.docker.com/desktop/, e.g. for [Linux](https://docs.docker.com/desktop/setup/install/linux/ubuntu/#install-docker-desktop)
  * Install Docker Engine (https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository)
    ```bash
    # Add Docker's official GPG key:
    sudo apt update
    sudo apt install ca-certificates curl
    sudo install -m 0755 -d /etc/apt/keyrings
    sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
    sudo chmod a+r /etc/apt/keyrings/docker.asc
    
    # Add the repository to Apt sources:
    sudo tee /etc/apt/sources.list.d/docker.sources <<EOF
    Types: deb
    URIs: https://download.docker.com/linux/ubuntu
    Suites: $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}")
    Components: stable
    Signed-By: /etc/apt/keyrings/docker.asc
    EOF
    
    sudo apt update

    sudo apt install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

    sudo docker run hello-world
    ```
  * Download DEB package "docker-desktop-amd64.deb"
  *  Install the DEB package
     ```bash
     sudo apt-get update
     sudo apt install (path to...)/docker-desktop-amd64.deb
     ```
  * Launch Docker desktop (accept "Docker Subscription Service Agreement")
* In case of error message "The path ... is not shared from the host and is not known to Docker."
  * Docker Desktop → Settings (or Preferences) → Resources → File Sharing. Add the problematic path (or its parent directory) to the shared list
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

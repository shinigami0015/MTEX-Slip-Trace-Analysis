# MTEX-Slip-Trace-Analysis

This repository hosts an automated workflow for crystallographic slip trace identification, Schmid factor mapping, and deformation compatibility analysis in dual-phase ($\alpha$-HCP / $\beta$-BCC) alloys using MATLAB and the MTEX toolbox.

---

## 🛠️ Core Features

* **Data Cleanup:** Automatic Euler-to-spatial alignment, scanning artifact purging, and grain boundary smoothing.
* **Trace Analysis:** Direct overlay of theoretical surface slip traces on Inverse Pole Figure (IPF) maps.
* **Stress Mapping:** High-throughput calculation of maximum resolved shear stress (Schmid factors).
* **Local Compatibility:** Boundary-by-boundary computation and mapping of the geometric Luster-Morris ($m'$) parameter.
* **Aggregated Statistics:** Automated phase statistics compiling active slip family distributions across all processed files.

---

## 📁 Repository Structure

* **`main.m`** – *Master Script:* Automates loop runs, scans input folders, exports plots, and graphs aggregate stats.
* **`EBSD_processing.m`** – *Data Preparation:* Handles file import, angular corrections, grain sizing filters, and spline smoothing.
* **`match_slip_label.m`** – *Helper:* Maps table string annotations to native 3-index (BCC) or 4-index (HCP) Miller indices.
* **`lustermorris.m`** – *Compatibility Module:* Computes grain pair orientation vectors and visualizes boundary compliance maps ($m'$).

---

## 🚀 Quick Start

1. Download all repository `.m` files into a single local workspace folder.
2. Ensure **MATLAB** and the **MTEX Toolbox** are installed and added to your active path.
3. Open **`main.m`** and update your data paths to point to your local directories:
   ```matlab
   inputDir  = 'D:\YourPath\To\Input CTF';
   outputBaseDir = 'D:\YourPath\To\Output files';
# Density-Driven Optimal Control (D2OC)
### Python and MATLAB Implementation of Optimal Transport–Based Multi-Agent Control: Time-Averaged Trajectory Matching for Specified Reference Distributions

### [![Hugging Face Space Demo](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Spaces-yellow)](https://huggingface.co/spaces/D2OC/D2OC-Demo)



<p align="center">
  <img src="figs/simulation.gif" width="550">
</p>

This repository provides the official MATLAB implementation of  
**Density-Driven Optimal Control (D2OC)** — a novel multi-agent control framework based on  
**Optimal Transport (OT)** and **Wasserstein distance** for  
**non-uniform area coverage**, **multi-agent/multi-robot coordination**, and  
**fully decentralized multi-agent control**.

This code accompanies the following publication:

📄 **Paper (IEEE Transactions on Systems, Man, and Cybernetics: Systems)**  
👉 DOI: https://doi.org/10.1109/TSMC.2025.3622075
👉 arXiv: https://arxiv.org/abs/2511.12756


---

## 🚀 Overview

D2OC solves decentralized multi-agent area coverage by:

- modeling target maps as **probability densities**,  
- guiding agents using **Wasserstein-distance–driven OT potentials**,  
- applying linearized quadrotor-inspired dynamics,  
- computing control inputs via a **finite-horizon KKT-based MPC**,  
- enabling decentralization through **local weight sharing**.

Included in this repository:

- 8-state quadrotor LTI dynamics  
- OT-based target computation  
- decentralized weight update logic  
- Gaussian-mixture density fields  
- simulation data and parameters  
- live trajectory visualization  

---

📁 Repository Structure
The source code is organized into Python and MATLAB implementations:

## [Python]
Located in the /Python folder.
Contains the single file for Python-based implementation of the D2OC framework for cross-platform research and integration.
Just run D2OC_main.py


## [MATLAB]
Located in the /MATLAB folder.
Main_D2OC.m : Main simulation script.
environment/DF.mat : Reference density maps.
param/param07.mat : Control parameters.
update_weight_R2.m : Decentralized weight update rule.
hamilton_optimal_control... : OT-based target computation logic.

## ▶ How to Run (MATLAB)
1. Clone or download this repository  
2. Use MATLAB R2020a or later  
3. Ensure `sim_data/Sim_rev60.mat` exists  
4. Open `Main_D2OC.m`  
5. Select a test:

```matlab
cnt_sim = 2;   % use 2, 3, or 4
```
6. Run the script  
7. Watch live visualization of UAV trajectories and density evolution  

---

## 🔥 Key Features

- **Optimal Transport control** with Wasserstein distance  
- **Non-uniform density tracking**  
- **Decentralized multi-agent coverage**  
- **Lagrangian-based OT point selection**  
- **Finite-horizon KKT/MPC formulation**  
- **UAV-ready MATLAB implementation**  
- **Scalable to many agents**

---

## 📘 State Dynamics

State vector:

```
x = [x, x_dot, theta, theta_dot, y, y_dot, phi, phi_dot]'
```

where:  
- x, y = positions  
- theta, phi = pitch/roll angles  
- derivatives = velocities  

---

## 📌 Applications

- Search & Rescue (SAR)  
- Environmental monitoring  
- Gas plume / wildfire mapping  
- Persistent surveillance  
- Agricultural field scanning  
- Exploration & inspection  

---

## 📊 Keywords

- density-driven optimal control  
- optimal transport control  
- wasserstein distance  
- multi-agent systems  
- decentralized control  
- distributed UAV control  
- coverage control  
- non-uniform area coverage  
- mpc  
- optimal control  
- uav robotics  
- density control  
- multi-robot coordination  

---

## 📖 Citation

If you use this code, please cite:

```
S. Seo and K. Lee, “Density-Driven Optimal Control for Efficient and Collaborative Multiagent Nonuniform Coverage,”
IEEE Transactions on Systems, Man, and Cybernetics: Systems, 2025.
DOI: https://doi.org/10.1109/TSMC.2025.3622075
```

(arXiv version: https://arxiv.org/abs/2511.12756)

---

## 📝 License
Released under the **MIT License**.

---

## 📧 Contact
Maintained by **Kooktae Lee, Ph.D.**  
Associate Professor, Mechanical Engineering, New Mexico Tech

# Density-Driven Optimal Control (D2OC)
### MATLAB Implementation for Optimal Transport–Based Multi-Agent Coverage Control

<p align="center">
  <img src="figs/d2oc.png" width="550">
</p>

This repository provides the official MATLAB implementation of  
**Density-Driven Optimal Control (D2OC)** — a novel multi-agent control framework based on  
**Optimal Transport (OT)** and **Wasserstein distance** for  
**non-uniform area coverage**, **UAV coordination**, and  
**fully decentralized multi-agent control**.

This code accompanies the following publication:

📄 **Paper (IEEE Transactions on Systems, Man, and Cybernetics: Systems)**  
👉 DOI: https://doi.org/10.1109/TSMC.2025.3622075
*(arXiv version coming soon)*

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

## 📁 Repository Structure

```
Main_D2OC.m                 → Main simulation script
environment/DF.mat          → Reference density maps
param/param07.mat           → Control parameters
sim_data/Sim_rev60.mat      → Simulation configurations
update_weight_R2.m          → Decentralized weight update rule
hamilton_optimal_control... → OT-based target computation
```

---

## ▶ How to Run

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
K. Lee, “Density-Driven Optimal Control (D2OC) for Multi-Agent Systems,”
IEEE Transactions on Systems, Man, and Cybernetics: Systems, 2025.
DOI: https://doi.org/10.1109/TSMC.2025.3622075
```

(arXiv version will be added soon.)

---

## 📝 License
Released under the **MIT License**.

---

## 📧 Contact
Maintained by **Kooktae Lee, Ph.D.**  
Associate Professor, Mechanical Engineering, New Mexico Tech

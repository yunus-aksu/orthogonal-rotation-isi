# orthogonal-rotation-isi
# Beyond Varimax: A Topological Stress Test for Orthogonal Rotation Criteria

This repository contains the official, fully reproducible R replication framework for the manuscript **"Beyond Varimax: Monte Carlo Stress Test for Orthogonal Rotation Criteria"**. 

This research challenges the universal reliance on the Varimax method as a default choice in factor rotation. By analyzing **449 million simulation runs** across 616 distinct data topologies (varying sample sizes, correlation structures, and variable counts) on High-Performance Computing (HPC) systems, we evaluate the numerical stability and variance redistribution behaviors of ten orthogonal rotation methods. 

We introduce the **Interpretability-Stability Index (ISI)**—a novel, system-level framework combining **Global Dominance (GD)** and **Selection Decisiveness (SD)**—to help researchers systematically navigate the empirical Pareto trade-off between structural simplicity and information concentration.
## 📂 Repository Structure

To run the replication pipeline, organize your local repository directory as follows:

orthogonal-rotation-isi/
├── reproduction.R                 # Main replication script (generates all tables & figures)
├── data/
│   ├── MASTER_scenario_freq.rds   # Raw simulation frequencies (from TRUBA HPC outputs)
│   └── MASTER_scenario_winner.rds # Raw simulation winner logs (from TRUBA HPC outputs)
└── output/                        # Automatically created directory for exported results
    ├── table_1.csv                # Replicated Table 1 (Loading Gap & SSI Statistics)
    ├── table_2.csv                # Replicated Table 2 (ISI Weighted Sensitivity Analysis)
    ├── table_3.csv                # Replicated Table 3 (Empirical validation on H-S & IPIP datasets)
    ├── figure_1.png               # Replicated Figure 1 (High-resolution Pareto Frontier in PNG)
    └── figure_1.pdf               # Replicated Figure 1 (Vectorized Pareto Frontier in PDF)
Note on Data Availability: If the raw .rds simulation files are not placed in the ./data/ folder, the script will automatically execute a fallback mechanism using the exact published parameters to output and export the tables, ensuring the script never breaks.

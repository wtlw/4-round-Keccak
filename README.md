# 4-round-Keccak
The code in this repository is related to the paper "Exploiting Output Bits and the $\chi$ Operation in MitM Preimage Attacks on Keccak."

The *paper code* folder contains the experimental code for the paper, all written in C++. Each subfolder corresponds to specific results from the paper: `512_noweak` contains code for the first experimental result, `512_weak` corresponds to the second, and `384_weak` to the third. To use this code, you need to install the Gurobi solver, which is available at [Gurobi's website](https://www.gurobi.com). The 3-round Keccak-512 experiment is the result and the experimental code described in the article.


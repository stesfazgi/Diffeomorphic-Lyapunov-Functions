# Learning Geometrically-Informed Lyapunov Functions with Deep Diffeomorphic RBF Networks

This repository contains an implementation of the paper [1],[2]:

[1] Tesfazgi, S., Sprandl, L., & Hirche, S. "Learning Diffeomorphic Lyapunov Functions from Data" in International Conference on Machine Learning (ICML) - Workshop on Geometry-grounded Representation Learning and Generative Modeling, 2024 [full paper](https://openreview.net/pdf?id=DUjZJe7wfF)
[2] Tesfazgi, S., Sprandl, L., & Hirche, S. "Learning Geometrically-Informed Lyapunov Functions with Deep Diffeomorphic RBFNs" in International Conference on Artificial Intelligence and Statistics (AISTATS), 2025 

---

## Requirements
- Matlab (tested on R2022b)
- LASA dataset (create a folder called LASA and include contents from https://bitbucket.org/khansari/lasahandwritingdataset/src/master/)
---
## Structure

Some code is duplicated in multiple subfolders. The idea behind that is that the following subfolders can be used as standalone packages
- diffeo-lyapunov/lya-limit           Code for learning Lyapunov functions for limit cycle systems
- diffeo-lyapunov/lya-multi-eq        Code for learning Lyapunov functions for multi-equilibria systems
- diffeo-lyapunov/lya-single-eq       Code for learning Lyapunov functions for single-equilibrium systems

---
## References
If you found this software useful for your research, consider citing us.
```
@inproceedings{Tesfazgi2025,
title={Learning Geometrically-Informed Lyapunov Functions with Deep Diffeomorphic {RBF} Networks},
author={Tesfazgi, Samuel and Sprandl, Leonhard and Hirche, Sandra},
booktitle={The 28th International Conference on Artificial Intelligence and Statistics},
year={2025},
url={https://openreview.net/forum?id=LysApRPxJ8}
}
```

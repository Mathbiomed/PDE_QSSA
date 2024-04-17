# PDE QSSA

Matlab code for simulating quasi-steady state approximations (QSSAs) of biochemical systems assuming a PDE system.
Related manuscript : Beyond Homogeneity: Assessing the Validity of the Michaelis–Menten Rate Law in Spatially Heterogeneous Environments, submitted to PLoS Computational Biology

## Code Description
1. Fig2
> The codes in this folder are designed to compare two models for approximating simple enzyme kinetics networks: the standard quasi-steady state approximation (sQSSA) and the total quasi-steady state approximation (tQSSA). Each code assumes different initial conditions.

2. Fig3

> fig3bc_.m and fig3de_.m: These codes are designed to compare the amount of product produced when the enzyme is localized in a small area vs. when it is not.
Similar to the codes in the Fig2 folder, these codes compare the sQSSA and tQSSA models.\\
> fig3i_.m: This is code to observe how inaccurate sQSSA becomes depending on the heterogeneity of the intial conditon.

2. Fig4
> Code to analyze the performance of QSSAs in a more complex model, i.e. a model involving two enzymes.
Here we can observe how well the QSSA models describe the zero-order ultrasensitivity of the full model.
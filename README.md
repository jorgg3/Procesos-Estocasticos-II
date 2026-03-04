# Stochastic Processes II: Computational Implementations

![Language](https://img.shields.io/badge/Language-R-blue.svg)
![Tool](https://img.shields.io/badge/Tool-RMarkdown-orange.svg)
![Library](https://img.shields.io/badge/Library-ggplot2-green.svg)

## Overview
This repository contains R and RMarkdown simulations developed for the Stochastic Processes II course. The project focuses on the computational implementation of various stochastic models, empirical verification of theoretical limit theorems, and the estimation of ruin probabilities in risk theory.

## Features & Implemented Models

The core of this repository is divided into three main computational exercises:

1. **Renewal Processes**
   * Algorithm to simulate sample paths of general renewal processes.
   * Empirical verification of the Elementary Renewal Theorem.
   * Step-function visualizations of the counting process $N_t$.

2. **Compound Poisson Processes**
   * Simulation of arrival times (Poisson process) and jump sizes.
   * Vectorized generation of interarrival times for computational efficiency.
   * Visualization of the accumulated jump trajectory $X_t$.

3. **Cramer-Lundberg Risk Model**
   * Implementation of the classical risk process: $R_t = u + ct - \sum Y_n$.
   * Monte Carlo simulation to estimate the probability of ruin $\psi(u)$ using exponential claims.
   * Extension of the model to a general renewal process with Gamma-distributed interarrival times and Exponentially-distributed claims.
   * Verification of the Net Profit Condition.

## Technologies Used
* **R:** Core programming language for statistical computing and Monte Carlo simulations.
* **RMarkdown:** Used for weaving together the theoretical mathematical explanations (via LaTeX) and the executable code.
* **ggplot2:** Utilized for creating high-quality, professional visualizations of the stochastic trajectories.

## Getting Started

To run these simulations locally:

1. Clone this repository:
   ```bash
   git clone [https://github.com/](https://github.com/)[jogg3]/[Procesos-Estocasticos-II].git

# Hamiltonian Learning, utils and demos

Repository for [Practical and Efficient Hamiltonian Learning](https://arxiv.org/abs/2201.00190).

## Overview

### Oracles

The oracle contained in *noisy_oracle.jl* mimics the actual evolution of the quantum system, by brutal-force state evolution under certain Hamiltonian. Specifically
- *oracle_f* extracts the second-order pauli error rates, as defined by equation (8) and (11)
- *oracle_s* gives the stage 2 measurements, defined by equation (15)

### Utils

*utils.jl* implements separately all the relevant procedures described in [the article](https://arxiv.org/abs/2201.00190), for instance, bins detection and peeling process (**Fig 1. (b)**) are wrapped in doPeel().

### hamLearning_xxx

These files glue and invoke all the subroutines together and demonstrate our algorithm under different settings (different Hamiltonian, different noise level, etc)

### data folder

The reconstruction results are stored in json format. Each json file contains
- a dictionary, representing the original Hamiltonian's parameters, i.e., the $s$ in main article
- a list of reconstructed parameters, containing
  - the reconstructed Hamiltonian
  - two numbers, counting the calls to the oracles in two subsequent stages.

### plot folder

Data analysis (with python) and results.



## Selected Results

- **Figure 3** in our article, demonstrating the supression of noise
    ![error](plot/Ising_n=6_varb_boxplot.png)
    ![ising_varb](plot/Ising_n=6_varb_violin.png)

- **Figure 4 (a)** in the article, random TFIM
    ![randomIsing_varn](plot/strictRandomIsing_n%3D1-7_violin.svg)

- **Figure 4 (d)** in the article, estimated Hamiltonian for the $\text{H}_4$ (8 qubits) molecule
    ![Molecules](plot/H4_new.svg)
<!-- - random Ising, under various noise level:

    ![ising1](plot/RandomIsing_n=4_varNoise_violin.svg)

- Ising, various noise level:

![ising2](plot/Ising_n=4_varNoise_violin.svg)

## Ising, scaling plot for different qubit number

random:
![ising3](plot/strictRandomIsing_n=1-7_violin.svg)

without randomness
![ising4](plot/strictIsing_n=1-8_violin.svg)

## Random Ising, various b

![var_b](plot/Ising_n=4_varb_violin.svg)

## Chemical Ham

LiH4, top 25 largest terms

![topTerms](plot/LiH4_top_20_terms_barplot.svg)

reconstruction with various b

![var_b](plot/Ising_LiH4_n=6_varb_violin.svg)

![var_b_plot](plot/LiH4_top_25_terms_scatterplot.svg)

## Random Ising n=4, increasing taylor expansion & fitting order $O(t^m)$

only even order is included, i.e. $m=6$ indicates fitting at order $t^0, t^2, t^4, t^6$

![var_m_plot](plot/RandomIsing_n=4_var_m_violin.svg)

## other bar plots

Ising, n=6
![Isingbar](plot/Ising_n=6_top_20_terms_barplot.svg)

## axis color changed

![error](plot/Ising_n=4_varNoise_violin_axisWithColor.svg)

![error](plot/LiH4_top_27_terms_redo_barplot.svg)

## box plot, demonstrating suppression of noise

![error](plot/Ising_n=4_varb_boxplot.svg)

## other chemical Hamiltonians

### Hchain, 4 atoms, 8 qubit, b=6, strengh * 10

top 60 terms

![error](plot/H4_top60_terms_barplot.svg)

### Hchain, 3 atoms, 6 qubit, b=4, 5, 6, 7, strengh * 10

top 20 terms

![error](plot/H3_top_20_terms_barplot.svg)

![error](plot/H3_n=6_varb_violin.svg)

### Hchain, 2 atoms, strengh * 10

full 14 terms

![error](plot/H2_full14_terms_barplot.svg)

![error](plot/H2_n=4_varb_violin.svg)

### revised figure, Hchain with 4 atoms

![error](plot/H4_new.svg) -->

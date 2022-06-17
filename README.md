# Hamiltonian Learning, utils and demos

repository for [Practical and Efficient Hamiltonian Learning](https://arxiv.org/abs/2201.00190)

## noise

random Ising, various noise level:

![ising1](plot/RandomIsing_n=4_varNoise_violin.svg)

Ising, various noise level:

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

![error](plot/H4_new.svg)

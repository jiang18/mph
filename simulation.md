---
title: Simulation
layout: template
filename: simulation.md
---

## Simulating phenotypes based on a list of GRMs
```
mph --simulate --num_phenotypes 100 --grm_list chr.grms.txt --heritability 0.5 --output sim_pheno
```
In the *i*th row of the GRM list file are GRM(*i*) and VC(*i*). MPH simulates total genetic values (**g**) by sampling **g** from N(**0**,**V**) in which **V** is equal to the sum of all GRM(*i*)\*VC(*i*). MPH further simulates phenotypes by adding an error term (**e**) to **g** based on heritability.

# EcoliPersistence
The program code and data files for generating Figures 3D-3F and 4 in [the article by Umetani et al](https://doi.org/10.1101/2021.10.28.466227).

The “Data” folder includes three folders with the names of growth phases from which *Escherichia coli* cells were sampled:
- Exponential
- EarlyStationary
- LateStationary

Each folder contains several XLS files for the raw data from the single-cell analysis of *E. coli* cultured in the microfluidic device. The data contain the information of cell size (area), mean fluorescent intensity, and cell lineage connections.

The program code (main.c) reads the XLS files placed directly under the “Data” folder. Hence, copy all of the XLS files from either of the three folders with the names of growth phases to the “Data” folder for analysis. 

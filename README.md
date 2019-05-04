# Design and Optimization of Phase Change Material Composites for Thermal Management

This is the data and code used for a class project titled "Design and Optimization of Phase Change Material Composites for Thermal Management". The goal of this project is to identify the optimal composition of the composite material heatsink that maximizes the heat absorbed by the heatsink system at 250, 500, and 1000 seconds after a heat pulse is applied. 

For this objective, two initial data sets were created using Latin Hypercube Sampling in conjuction with a finite element analysis model representing
the physical system. The two data sets contain 500 points (preliminary data set) and 10,000 points (primary data set). The preliminary data set
is used to identify the relationship between different radially varying volume fractions of materials in the system. The primary data set is then 
used to determien the number of samples required to fit the optimizers being investigated. 

This response surface is used in our materials informatics approach to identify the optimal configuration of the materials system while minimizing
the number of experiments (or calculations) requried. This surface will also be useful in generating a forward approach to future heatsink design
with this materials. 

This repository contains the two original datasets mentioned above, as well as, Matlab
codes used to for the optimization. 

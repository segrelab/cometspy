# COMETS Capabilities

The core of COMETS is the simulation of the growth of microbial populations by maximizing an objective reaction, usually biomass production, in an iterative way. An initial amount of biomass of one or more species (defined by their metaboli models) is seeded in an environment containing a list of specified nutrients. In each iteration, both the amount of biomass and the environment are updated using FBA predictions. 


## Capabilities in space

COMETS is capable of simulating microbial growth in a spatially structured environment. This is achieved by partitioning the simulation space ("world") in a "grid" of smaller spaces. Inside each space, growth is considered to be well-mixed. Both nutrients and biomass can then propagate to contiguous spaces as simulation proceeds. 

- 2D and 3D simulation "worlds"

In addition to well-mixed conditions, COMETS can simulate 2D and 3D spatially structured environments. This enables simulation of, for instance, growth of colonies on 2D surfaces such as a Petri dish, or 3D structures such as tumors, bacterial colonies in 3D matrix, etc.


- Diffusive and convective propagation of biomass in space

In simulations with spatial structure, different modes of biomass propagation are implemented. The diffusive mode simulates the propagation of free swimming motile bacteria, while the convective mode simulates the propagation of bacteria by mutual pushing. The two modes of propagation can be combined.


- Substrate-dependent nutrient and biomass propagation

The diffusivity of nutrients as well as the propagation properties of the biomass depend on the substrate properties such as agar density, cell substrate friction coefficient etc. 


- Boundary conditions

Two types of boundary conditions are implemented. Fixed value and fixed source or sink rate. 


## Biological capabilities

COMETS features many interesting biological capabilities to refine and improve the predictions of stoichiometric models, as well as to simulate different types of biologically realistic conditions. 

- Lag-phases in microbial growth

Lag-phases are modelled as simulated growth activation of the colonies.


- Continuous (chemostat) and batch growth modes

In chemostat mode, the user controls the rate of replenishment of the nutrient. In batch mode, the user controls dilution and frequency.


* Simulation of multispecies communities

Simulation of two or more species (up to hundreds) can be performed in both species overlapping or non-overlapping spatial distribution, or in well-mixed conditions. 


* Parsimonious dFBA

Usually, any metabolic model has multiple optimal solutions. One way to choose among them is to assume that the cell will minimize the total flux through the metabolic network. To achieve this, parsimonious FBA first optimizes the objective function, e.g. growth. Then, it performs a second optimization by fixing growth at the previously obtained optimal level and minimizing the total flux through the network. 


* Cell death

A simple model of cell death is implemented with each species assigned death rate.


- Neutral population drift

The presence of demographic noise can result in random variations in the abundance of different species in a simulation. This is especially useful in the batch-growth mode, where dilution bottlenecks can have a significant impact on growing populations. 


* Evolutionary processes

Comets allows for evolutionary simulations, including random mutation and drift during simulations. Currently, the only mutations that are available are reaction deletions. 


## Computational capabilities
COMETS software is implemented in the JAVA language. Therefore, it is highly portable and independent on the operative system. COMETS offers the following simulation capabilities. 

- Graphical User Interface (GUI) 

In addition to the command line, COMETS simulations can be run using a graphical user interface that includes visualization tools. 


- Parallelized dFBA 

Runs in multi-CPU systems as multi-threaded process for greater computational performance.


- MATLAB toolbox

A toolbox in MATLAB for modifying the input files for COMETS in a programmatic way. 


- Python toolbox

A toolbox in Python for modifying the input files for COMETS in a programmatic way.

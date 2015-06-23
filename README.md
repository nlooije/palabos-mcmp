# palabos-mcmp
palabos multi-component, multi-phase template codes

# extensions to the palabos code base:
The include directory includes files for:
* arbitrary equations of state using the Yuan and Schaeffer model
* different force implementation (Guo, SC, EDM) through a generalized model by Lycett-Brown
* modified data processors for simulating rising bubbles in a periodic domain

# how to use
Download the ./include folder and set its path in the Makefile. Include the relevant files and use the new function calls in the simulation case and compile/run as normal

# tkeBudget
OpenFOAM's functionObject used to calculate turbulence kinetic energy budget of incompressible flows 

Turbulence kinetic energy equation given as:

$ [\frac{\partial}{\partial t}+U_j\frac{\partial}{\partial x_j}]k = \frac{\partial}{\partial x_j}\{-\frac{1}{\rho}\langle pu_i\rangle\delta_{ij}-\frac{1}{2}\langle q^2u_j\rangle+2\nu\langle s_{ij}u_i\rangle\}-\langle u_iu_j\rangle \frac{\partial U_i}{x_j}-2\nu\langle s_{ij}s_{ij}\rangle $

where:
$U$ denotes mean velovity; 
$k=\frac{1}{2}\langle q^2\rangle=\frac{1}{2}\langle u_iu_i\rangle$ denotes average fluctuating kinetic energy per unit mass;
$s_{ij} = \frac{1}{2}(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i})$ denotes symmetric part of the strain rate tensor;

Therefore,each term of tkeBudget can be evaluate.
For a detailed introduction, refer to [Introduction to turbulence/Turbulence kinetic energy](https://www.cfd-online.com/Wiki/Introduction_to_turbulence/Turbulence_kinetic_energy).

# Introduction

This project deals with the thermo-elastic and limit analysis of reinforced concrete thin plates. Extension to other materials are possible by changing the material properties in the study. Here is the [reference](https://hal.archives-ouvertes.fr/tel-02297575/).

# Temperature gradient

The temperature gradient across the thickness of plate is calculated by solving the one dimensional heat diffusion problem.

Thermal properties of concrete, such as The heat capacity, conduction, water content, are obtained according to the Eurocode. Temperature distribution are calculated or different plate thickness, from 10cm to 30cm, and stored in files at every 30mins.

# Material properties

## Elastic properties
Only the elastic properties of concrete will be considered due to the fact that the volumetric fracture of steel is small. The degradation of elastic properties of concrete is obtained according to the EuroCode.
## Strength properties
Both concrete and steel bars are accounted in the calculating the strength domain of reinforced concrete plate sections. The degradation of the concrete and steel properties at elevated temperature is obtained by the EuroCode.
# Thermo-elastic deformation


## Kirchhoff-Love plate

Analytical solution for thermo-elastic deformation of Kirchhoff-Love plate
The temperature gradient is uniform across the plate surface.

Three different cases with different geometry and boundary conditions:
- Case 1: Rectangular plates simply supported on four sides.
- Case 2: Rectangular plates simply supported on two opposite sides and free on the other two opposite sides.
- Case 3: Rectangular plates with uniformly distributed joints

## Von Karman plate

Semi-analytical solution based on the Galerkin method for rectangular walls subjected to both the self-weight and a temperature gradient across the thickness of walls.

Two different cases with different geometry and boundary conditions:

- Case 1: Rectangular plates simply supported on four sides.
- Case 2: Rectangular plates simply supported on two opposite sides and free on the other two opposite sides.

# Strength domain of reinforced section
## Lower bound by a static approach

## Upper bound by a kinematic approach

# Limit analysis

## Lower bound by a static approach 

## Upper bound by a kinematic approach 
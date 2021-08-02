# Scaling tests with the graphene + proton example

## Getting started

I have written a script that generates an input file for an SCF calculation of the ground state of a tiling of the graphene lattice.
You can run this script as follows:

python3 generate_graphene_proton_supercell.py 2 2

will generate the 2x2 supercell that you've been studying.
By substituting the two arguments (here, 2 and 2) you can change the shape of the tiling.
E.g., 4 3 would tile the cell 4 times along the first lattice vector and 3 times along the second lattice vector.
The total number of atoms will be 1 (proton) plus 2 times the product of the two arguments.
E.g., 2 2 will give you 1 + (2 atoms)*(2 cells in direction A)*(2 cells in direction B).

## To do list

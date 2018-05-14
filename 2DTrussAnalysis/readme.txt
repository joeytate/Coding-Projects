This code takes the given coordinates of nodes on a truss, the area and Young's Modulus 
of elements between each connected node, any the location, direction, and magnitude of any applied forces
on the nodes, and the known displacements of the nodes to calculate the missing reaction forces and displacements
of the linear system using the Penalty Method. It then plots the resulting truss with its proper displacements.

project_input.txt can be altered to change the number of nodes and elements on the truss, the
size and Young's Modulus of the elements, as well as force and displacement values at each node.
# cst-element-fortran
CST elements analyzing using finite element method in FORTRAN

== Analyzing CST elements using FORTRAN90 ==

==== Developed by Mohammad Sajednia ====

Inputs:

1. nodes.txt
First row: the number of the nodes
Column 1: X coordinate of the node
Column 2: Y coordinate of the node

2. elements.txt
First row: the number of the elements
Column 1: the number of the first node
Column 2: the number of the second node
Column 3: the number of the third node

3. boundaries.txt
First row: the number of the constraints
Column 1: the number of the constrained node
Column 2: the direction of the constraint (1 for X, 2 for Y)
Column 3: the value of the constraint (0 if the constraint is a support and any value for support settlements)

4. nloads.txt
First row: the number of the nodal loads
Column 1: the number of the node
Column 2: the magnitude of the nodal force
Column 3: the angle between the force vector and the positive direction of the X axis

5. tractions.txt
First row: the number of the tractions
Column 1: the number of the first node
Column 2: the number of the second node
Column 3: the magnitude of the traction between two nodes (on the edge)
Column 4: the angle between the traction vector and the positive direction of the X axis

6. vloads.txt
First row: the number of the body forces
Column 1: the number of the element that includes the body force
Column 2: the magnitude of the body force
Column 3: the angle between the body force vector and the positive direction of the X axis

Note 1: Don't define negative angles (Ex: define 270 degree instead of -90 degree)

Note 2: When you run the program, it'll ask you if there is at least one nodal force, traction, or body force. Type "1" and ENTER if the answer is YES.

Note 3: After running the program, a file named "Outputs.txt" in the same directory and the results will be printed on it.

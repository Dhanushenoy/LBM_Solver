# LBM_Solver
Lattice Boltzmann method based solver for incompressible Navier-Stokes equation
*Currently 2d with D2Q9 method, no turbulent model yet
*currently solving just the velocity field

In linux, just execute the "run" file. 
The boundary conditions, collision step, stream step etc are in the folder SRC
Currently, the generated .o and .mod files are moved everytime after creating (Check run file)
VTK file can be opened with Paraview, or there is also a subroutine for simple x,y,u,v output

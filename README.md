# The Heat Method for Distance Computation
I have implemented the paper "The Heat Method for Distance Computation" on triangle meshes. For further details it is worth visiting 
[Keenan Crane's Page](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/).


## Paper summary

The paper proposes the heat method for solving the single- or multiple-source shortest path problem on both flat and curved domains. 

It splits the distance computation procedure into two stages: 
- First find the direction along which distance is increasing by integrating the heat flow for a fixed time and evaluating the vector field.
- Compute the distance itself by solving a Poisson equation.

The method is based on solving a pair of sparse linear systems. Hence, it is easy to implement. Moreover, the sparse systems can be factorized once and
reused. 


## Some Implementation Details

I have used [Geometry Central](https://geometry-central.net/) and [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) libraries for their data structures, operators, and solvers. For visualization,
I have used [Polyscope](https://polyscope.run/)



## Benchmarks
| **Mesh**   | **Number of Triangles** | **Setup Time**        | **Solving Time**       | 
|---------------|--------------|------------------------|-------------------------------|
| Cat           | $1k$         | $0.003 s$             | $0.0003s$                      | 
| Man           | $24k$        | $0.05 s$              | $0.01s$                        |
| Centaur       | $31k$        | $0.06 s$              | $0.01s$                        |
| Dragon (Simplified)| $64k$   | $0.2 s$               | $0.04s$                        |
| Bunny         | $69k$        | $0.331 s$             | $0.04s$                        |
| Lucy          | $100k$       | $0.336 s$             | $0.05s$                        |
| Strawberry    | $314k$       | $1.67 s$              | $0.13s$                     |

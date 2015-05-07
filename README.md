# gradientgenmodeling

Fortran code used to numerically solve the 3D concentration profile of microfluidic gradient generators used in the manuscript:
  
Gorman, B.R. & Wikswo, J.P. Characterization of transport in microfluidic gradient generators. Microfluid Nanofluid 4, 273-285 (2008). [doi: 10.1007/s10404-007-0169-0] (http://dx.doi.org/10.1007/s10404-007-0169-0).  
  
It is an implementation of the line-by-line method discussed in Suhas Patankar's book, [Numerical Heat Transfer and Fluid Flow](http://www.ewp.rpi.edu/hartford/~ernesto/F2012/CFD/Readings/Patankar-NHTFF-1980.pdf) (1980), page 64-66. Essentially, it is a combination of Gauss-Seidel and the Thomas algorithm.  In other words, the algorithm solves directly along a line using the Thomas algorithm while holding the other variables explicit, and then sweeps the line along the domain.  Then, it switches the direction - i.e., switches from sweeping left-to-right to sweeping top-to-bottom. Thus, it is an efficient way to bring information from the peripheral parts of the domain to the center and vice-versa.

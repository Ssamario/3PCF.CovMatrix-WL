# 3PCF.CovMatrix-WL

Code to compute the Gaussian piece of the 3PCF in a harmonic basis of scalar fields over the sphere within Limber approximation
Authors:
Sofia Samario (ssamario@icf.unam.mx)
Alejandro Aviles (avilescervantes@gmail.com)

To perform the integration routines, we take advantage of the GNU Scientific Library (GSL) package that provides us with the Bessel functions.

To run the code for the first time 
Run in the terminal 
gcc -Wall -I/path/to/gsl/include -c IntegralCovMatrix.c
then 
gcc -Wall -I/path/to/gsl/lib IntegralCovMatrix.o -lgsl -lgslcblas -lm






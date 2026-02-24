# 3PCF.CovMatrix-WL

Code to compute the Gaussian piece of the 3PCF in a harmonic basis of scalar fields over the sphere within Limber approximation
Authors:
Sofia Samario (ssamario@icf.unam.mx)
Alejandro Aviles (avilescervantes@gmail.com)

This code uses the **GNU Scientific Library (GSL)** for the numerical integration routines and the evaluation of Bessel functions.

Make sure GSL is installed on your system before compiling.

### Compilation

To compile the code for the first time, run:
```
gcc -Wall -I/path/to/gsl/include -c IntegralCovMatrix.c
gcc -Wall IntegralCovMatrix.o -L/path/to/gsl/lib -lgsl -lgslcblas -lm -o IntegralCovMatrix
```

Replace ```/path/to/gsl/``` with the actual installation path of GSL on your system.

### Running the Code

After compilation, execute:
```
./IntegralCovMatrix
```

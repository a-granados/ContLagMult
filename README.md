<h4> Continuation of curves using Lagrange Multipliers </h4>
This software computes a curve defined implicitly by a function F from R^m to R^n, with n=m+1. That is, we have n equations and n+1 unknowns. Assuming no degeneracy, this deifnes a curves in R^n.
Soon more details about the theoretical background

You will need GSL libraries installed.

How to compile:
g++ -o execname cont_LagrangMult.c rk78.c Fun.c lu.c memoria.c -lgsl -lgslcblas

Fun.c is a file that you need to create. There you define the function F and its Jacobian. See the folder "examples" for examples.

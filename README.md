# ligfit

A program for fitting equilibrium binding data for homodimerization around
homo-bifunctional ligands.

The mathematical foundation for this project is laid out in the publication by
Mack et al.:

Mack ET, Perez-Castillejos R, Suo Z, Whitesides GM. Exact Analysis of
Ligand-Induced Dimerization of Monomeric Receptors. Analytical chemistry. 
2008;80(14):5550-5555. doi:10.1021/ac800578w.

This code originally written by Paul Barton while a PhD. candidate in Biological
Sciences at Carnegie Mellon University.

The seminal paper for fluorogen-activating single chain antibody fragments as
fluorogen-activation proteins: doi:10.1038/nbt1368
Studying the kinetics of binding for these systems provided the original
motivation for creating this project.

## An Introduction

This binding model possesses the distinct advantage of having an analytical
solution, providing the ability to exactly (without simplifying approximations)
compute how the system should behave theoretically under varied conditions.
However the computation of this model faces a particular challenge (one readily
demonstrated in the notebook of this project) in that it is highly sensitive to
floating point imprecision depending on the parameters of the system/experiment.
Put another way, with the right circumstances, a computation of the model may
yield a wildly inaccurate answer or fail utterly.

The response to this challenge involves utilizing ever more precise
representations of floating point numbers such that the imprecision becomes
irrelevant to the system, and its parameters, at hand. With this tactic in mind
I began to write the code to fit experimental data to this binding model using
SciPy and mpmath. The first providing an tools for solving error minimization
problems in curve fitting, and the latter allowing arbitrary levels of precision
in floating point numerical operations. This ultimately failed because, as I
came to learn, SciPy was only providing a wrapper around various mathematical
packages, often compiled fortran, which were not readily able to interact with
mpmath.

It then appeared that a reasonable strategy would be to adapt the fortran
packages I would need to employ for direct use with arbitrary floating point
precision. My research into this problem led me to employ MPFUN2015 by David H.
Bailey (http://www.davidhbailey.com/dhbsoftware/) as I adapted the necessary
subroutines from MINPACK so that they might successfully fit experimental data
with the model.

## Overall Project Composition and Requirements

Fortran lies at the heart of this project and I am using Fortran90/95. MPFUN2015
is a Fortran90 package, and the MINPACK from which I worked is also Fortran90.
To successfully use this project, you will need an adequate Fortran compiler.

Python plays an important role in this project, I am using Python3 and a myriad
of packages not in the standard library for various jobs. To use the core
function of ligfit (fitting data!), you will need NumPy and docopt. I won't
remark on the dependencies of NumPy here, but I will note that as of the time
of writing (2015-05-24) I had to install the traditional way rather than using
pip, as installing via pip did not also install the vital tool f2py. The
traditional way being: first download the source, then compile/install with
"python(3) setup.py install".

To make use of the notebook and the other various adjuncts I have provided, you
will need to install ipython, scipy, and matplotlib. Pip for these should do,
the setup.py installation script for this project should automatically get these
dependencies for you.



The Radial Basis Function Toolbox (RBFT) is a collection of functions for implementing RBF interpolation methods and RBF methods for the numerical solution of PDEs on scattered centers located in complexly shaped domains.   The toolbox is available in for Matlab (MRBFT) and a Python (PRBFT) version will be released in the future.  The Matlab version uses the Multiprecision Computing Toolbox (http://www.advanpix.com/) to seamlessly implement extended precision floating point arithmetic in all RBFT routines. 

Comments, questions, bug reports, code requests, etc. can be sent to sarra@marshall.edu


The functionality of the toolbox is organized via object oriented programming into several classes:

    rbfX - basic RBF method functionality
        gax - Gaussian RBF
        iqx - Inverse Quadratic RBF
    rbfCenters - center locations
    rbfCentro - reduced flop count and storage algorithms for RBF methods in symmetric domains.
    Functions - test functions and derivatives.

The toolbox comes with a collection of scripts that demonstrate its usage, benchmark its performance, and verify that its algorithms produce the correct results.  The scripts are located in the following folders.

    \examples
    \tests
    \benchmarks


If the RBFT has been significant to a project that leads to an academic publication, please acknowledge that fact by citing the project.  The academic reference for the RBFT is this paper (http://www.scottsarra.org/math/papers/sarraMRBFT.pdf).  The BibTex entry for the paper is

@Article{Sarra2016,
  Author    = {S. A. Sarra},
  Title     = {The {M}atlab Radial Basis Function Toolkit},
  Journal   = {Journal of Open Research Software, under review},
  year      = 2016,
  url       = "www.scottsarra.org/rbf/rbf.html",
}

or in plain text:

S. A. Sarra.  The Matlab Radial Basis Function Toolkit.  Under review, Journal of Open Research Software, 2016. 

Thank you!



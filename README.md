Irreducible angular momentum and spin eigenspaces on atomic subshells
=====================================================================

Mathematica implementation of the algorithm described in the paper **Efficient algorithm for many-electron angular momentum and spin diagonalization on atomic subshells**, see Ref. 1 below.

The main code resides in the Mathematica package `irredLSbase.m`, and the computations are performed in the notebook `irredLS.nb`. This notebook requires the **FermiFab** toolbox (see Ref. 2) and stores the results in `irredLS.m`.

The files `fermi2latex.m` and `irredLSlatex.nb` generate the Latex tables of the eigenspaces in the *tables* subfolder.

License
-------
Copyright (c) 2009-2014, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License
http://www.opensource.org/licenses/bsd-license.php


References
----------
1. Christian B. Mendl  
   Efficient algorithm for many-electron angular momentum and spin diagonalization on atomic subshells  
   preprint [arXiv:1409.6860](http://arxiv.org/abs/1409.6860) (2014)
2. Christian B. Mendl  
   The FermiFab toolbox for fermionic many-particle quantum systems  
   Comput. Phys. Commun. 182, 1327-1337 (2011), [arXiv:1103.0872](http://arxiv.org/abs/1103.0872)  
   URL http://sourceforge.net/projects/fermifab

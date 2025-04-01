About
=====

Associated with the publication of `Smoothed Particle Hydrodynamics and Magnetohydrodynamics <http://adsabs.harvard.edu/abs/2012JCoPh.231..759P>`_ in J. Comp. Phys., I made my "development" SPH code, NDSPMHD, public.

Features include:

- Full implementations of 1D, 2D and 3D hydrodynamics and magnetohydrodynamics as described in Price (2012)
- (v2.1) Implementations of two fluid and one fluid dust-gas algorithms and test problems, as described in Laibe & Price (2012a,b,2014a,b) and Price & Laibe (2015)

Note, however, that ndspmhd is *not* meant as a "production" SPH code in 3D, since much better codes exist
for this purpose (e.g. my own `PHANTOM <https://github.com/danieljprice/phantom>`_ code). ndspmhd is not parallel nor particularly optimised and is meant as a code for algorithmic experimentation, not production runs. 

Licence/conditions of use
--------------------------

The code is distributed under the GNU general public licence (v2.0). The only conditions of use aside from this are that you:

- Cite the paper (Price D.J, 2012, J. Comp. Phys. 231, 759-794, `arXiv:1012.1885 <http://www.arxiv.org/abs/1012.1885>`_) if you publish anything based on the code;
- Cite the `dust papers <references>`_ if you use the dust algorithms (and likewise for MHD); and
- Kindly send me a copy of any such manuscript *prior to acceptance* (i.e., on submission to a journal/proceedings). 

History
--------

ndspmhd was developed as part of my PhD research at the University of Cambridge from 2001-2004. It was used in nearly every paper I published during a 10-year period afterwards and contains complete working versions of most of the algorithms described in print in those papers.

The history of the public code is as follows:

- 08/05/2015: v2.1: Minor bug fix with build (thanks to Victor Moral)
- 05/05/2015: v2.1: Update including dust diffusion and non-ideal MHD algorithms and test problems.
- 26/06/2014: v2.0: Bug fix with release tarball (thanks to Marc Joos).
- 21/02/2014: v2.0: Major update including dust algorithms.
- 21/12/2010: v1.0.1: Minor bug fix with build.
- 26/10/2010: v1.0: First public version of ndspmhd posted on web (only minor changes from ASTROSIM version).
- 07/09/2010: A (non-GPL'd) version of ndspmhd was posted on my web page for the ASTROSIM summer school. 

Support/feedback
-----------------

Bug reports and feedback are always appreciated. 
.. _examples:

Running the examples
=====================

Preliminaries
-------------

After unpacking the tar file:

.. code-block:: bash

   tar xvfz ndspmhd-examples.tar.gz

and entering the examples directory,

.. code-block:: bash

   cd ndspmhd-examples/

You will see three subdirectories for hydrodynamical, MHD (both ideal and non-ideal) and dusty-gas tests:

.. code-block:: bash

   $ ls
   DUST/  HYDRO/  MHD/  MHD-NONIDEAL/

For hydrodynamics (in the HYDRO directory), there are 7 tests:

.. code-block:: bash

   test1-2Drandom-settling/  test3.1-2Dsodshock/  test6-2Dkhinstability/
   test2-2Dpairing-instability/  test4-1Drarefaction/  test7-3Dsedov-blast-wave/
   test3-1Dsodshock/  test5-1Dblastwave/

Whilst for MHD there are 7:

.. code-block:: bash

   test1-1.5Dbriowu-shock/  test3-1.5Dtoth-strong-mhdshock/  test5-2Dmhdrotor/  test7-2Dcurrent-loop-advection/
   test2-1.75Dmhd-shock/  test4-2.5Dalfvenwave/  test6-2Dorszag-tang-vortex/

and for dust-gas (in the DUST directory) there are 5:

.. code-block:: bash

   test1_dustybox/  test2_dustywave/  test3_dustyshock/  test4_dustdiffusion/  test5_dustsettle/

These are all standard test problems in the computational astrophysics literature and appropriate
references can be found alongside these tests in my `papers <https://ui.adsabs.harvard.edu/public-libraries/bRgpY8A3SYGxgGd3LOMQAw>`_ on SPH and SPMHD algorithms. The :ref:`reference list <references>` gives the most useful of these.

Running the tests
-----------------

It is preferable to go through the examples in the order given. For each test there is a README-blah.txt file in the subdirectory describing how to run the code step-by-step, so things are fairly self-explanatory.

For example, entering the first test directory

.. code-block:: bash

   cd test1-2Drandom-settling/

you will find the following files:

.. code-block:: bash

   $ ls
   Makefile  random.in  rpsph.in  splash.limits
   README-settling.txt  randomav.in  splash.defaults

...so simply follow the instructions in the README-settling.txt file...

.. code-block:: bash

   $ more README-settling.txt
   This is a test showing how a random particle distribution will relax to an
   ordered arrangement because of the Hamiltonian nature of SPH.

   Type "make clean" and "make" to re-compile the 2DSPMHD binary.

   Run the code using ./2DSPMHD random.in
   Or with output to a file: ./2DSPMHD random.in >& random.out &

   Plot the results using "nsplash random_0*.dat"
   ...

Note that you will need to have `SPLASH <https://splash-viz.readthedocs.io>`_ installed to view the code output.

Next steps...
-------------

After running a few of the examples, you will have a rough idea of how the code works and what the important input options are. From there you can proceed to explore the tests further or continue towards :ref:`setting up your own problem <tutorial>`... 
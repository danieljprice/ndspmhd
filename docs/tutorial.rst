.. _tutorial:

Running your own setups
=======================

Creating a new run directory
----------------------------

Compiling the code will produce a binary called either 1DSPMHD, 2DSPMHD or 3DSPMHD in the root-level ndspmhd directory. The procedure I use for starting a calculation is generally as follows:

Create a new directory and enter it:

.. code-block:: bash

   mkdir test-run
   cd test-run

Use the script in the ndspmhd/scripts directory to write a local Makefile:

.. code-block:: bash

   ~/where-I-put-it/ndspmhd/scripts/writemake.sh 1D > Makefile

then, after setting the location of the code via the NDSPMHD_DIR environment variable (best to set this permanently by putting it in your .bash_profile or equivalent)

.. code-block:: bash

   export NDSPMHD_DIR=~/where-I-put-it/ndspmhd

Now you can simply type

.. code-block:: bash

   make

which should compile the code and copy the binary to the current directory. A directory listing shows the binary and the local Makefile:

.. code-block:: bash

   $ ls
   1DSPMHD*  Makefile

Creating an input file
----------------------

Create an input file with default options by running the code with the name of an input file as the argument:

.. code-block:: bash

   $ ./1DSPMHD test.in
    number of runs =  1
    run =  1  runname = test.in
    input file test.in              not found
     would you like to create one with default options?
   y
    input file test.in              created successfully
   STOP exiting...

Editing the particle setup / initial conditions
-----------------------------------------------

When run with an input file as the argument, the code runs from scratch with the particle setup specified by the SETUP1D, SETUP2D or SETUP3D variable in src/Makefile. To edit the currently compiled setup,  use

.. code-block:: bash

   make edit

which will load the relevant file into your favourite editor (requires setting the EDITOR environment variable to be your favourite text editor). Alternatively you can make a new setup_blah.f90 file and specify this using SETUP1D=setup_blah.f90 in the Makefile. Then

.. code-block:: bash

   make

to rebuild the code in the local directory with the edited file.

You can also specify the setup in the local Makefile by giving an extra argument to the writemake script, i.e.,

.. code-block:: bash

   ~/where-I-put-it/ndspmhd/scripts/writemake.sh 1D setup_blah.f90 > Makefile

Then when you type "make" the code will always build using setup_blah.f90 as the initial conditions. This is how the Makefiles in the `examples <examples>`_ directories have been written.

Running parameter sweeps
------------------------

With more than one input file on the command line, e.g.

.. code-block:: bash

   ./1DSPMHD test1.in test2.in test3.in

the code will run each test concurrently (provided the previous test finishes normally). This is a convenient way of performing parameter sweeps (in serial) without any scripting.

Restarting the code
-------------------

It is possible to restart the code from the position of any binary dump file by simply giving the name of the file instead of the input file on the command line, e.g.

.. code-block:: bash

   ./1DSPMHD test_00020.dat

which will look for input options in a file called test.in. 
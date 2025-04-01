Installing
===========

System requirements
-------------------

The basic requirements for compiling and running ndspmhd are simply a reasonably modern Fortran compiler. The instructions below assume a unix-based operating system (e.g., Linux or Mac OS/X). There is no reason in principle why ndspmhd should not run under Windows, but I can offer no support in this regard (likewise for SPLASH).

Getting a Fortran compiler
---------------------------

As NDSPMHD is written in Fortran 90, you will need to have a modern Fortran compiler installed. Options include:

- `gfortran <http://gcc.gnu.org/wiki/GFortran>`_, the free Gnu Fortran Compiler, download binary version `here <http://gcc.gnu.org/wiki/GFortranBinaries>`_.
- `ifort <http://software.intel.com/en-us/intel-compilers/>`_, one of the most widely available commercial compilers (and is very good).

Download the code file
----------------------
.. list-table::
   :widths: 25 50 25
   :header-rows: 1

   * - File
     - Description
     - Date
   * - `ndspmhd-v2.1.tar.gz <https://users.monash.edu.au/~dprice/ndspmhd/ndspmhd-v2.1.tar.gz>`_ (143Kb)
     - The code
     - 08/05/2015
   * - `ndspmhd-examples-v2.1.tar.gz <https://users.monash.edu.au/~dprice/ndspmhd/ndspmhd-examples-v2.1.tar.gz>`_ (184Kb)
     - Exercises and examples associated with the code.
     - 05/05/2015

You will also want to download and install `SPLASH <https://github.com/danieljprice/splash>`_ to be able to look at the binary output files (in v3.x of splash, use "splash -f ndspmhd" to read this format; in v2.x of splash, "nsplash" reads this format). To correctly read and visualise the dust examples in ndspmhd v2.1 you will need v2.5.2 of SPLASH or later. 

Compiling the code
------------------

Uncompress the tar file:

.. code-block:: bash

   tar xvfz ndspmhd-v2.0.tar.gz

and enter the directory

.. code-block:: bash

   cd ndspmhd

Compile the code by typing "make" whilst specifying the number of spatial dimensions, e.g.

.. code-block:: bash

   make 1D

which gives the following error

.. code-block:: text

   make: WARNING: value of SYSTEM =  not recognised...
   =>set the environment variable SYSTEM to one listed
     in the Makefile and try again

   I suggest one of the following, based on detected Fortran compilers...

   make SYSTEM=gfortran
   (end of possible selections)

settings for different SYSTEMs are listed in src/Makefile. Either choose one suggested above or add your own customised setup in src/Makefile. Compile for your chosen system using (e.g.)

.. code-block:: bash

   make SYSTEM=gfortran 2D

or by setting SYSTEM as an environment variable, e.g. in bash:

.. code-block:: bash

   export SYSTEM=gfortran;
   make 2D

Compiling the code will produce a binary called either 1DSPMHD, 2DSPMHD or 3DSPMHD in the root-level ndspmhd directory.

Running the exercises...
------------------------

Once the code has compiled successfully, you should proceed to :ref:`running the examples <examples>`. 
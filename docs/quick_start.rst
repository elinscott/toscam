Quick start
===========

Performing a DFT+DMFT calculation with ``ONETEP+TOSCAM`` involves

1. Running a standard DFT singlepoint calculation using ``ONETEP``
2. Running DMFT on the DFT output using ``TOSCAM``

This guide will assume the reader is comfortable with performing the first step, and we will focus on the second. (For guidance on using ``ONETEP``, refer to the `ONETEP website <https://onetep.org>`_)

Required files
--------------
The following files must be contained in the working directory

input_onetep_dmft.txt
   This file specifies the settings for the DMFT calculation. Each setting is listed on a row of the file, in the format ``keyword=value``. Comments are to be prefixed by ``!``. Keywords are listed in :ref:`Page Keywords`

<seed>.dat
   The standard ONETEP ``.dat`` file. Note that this may be overwritten during the DMFT calculation; if a setting is specified both here and in ``input_onetep_dmft.txt`` the latter's contents take precedence.

<seed>.dkn
   Density kernel file from the initial DFT calculation

<seed>.tightbox_ngwfs
   NGWFs from the initial DFT calculation

mask_projections
   A file specifying the Hubbard projections to use as indicated by a string of ``T``'s and ``F``'s. For example, if one  wants to explicitly treat all five `3d` orbitals of a single-Hubbard-atom system, this file should read ``TTTTT``

Parallelism
-----------
This code has hybrid MPI-OpenMP parallelism. However, the number of MPI tasks for the DMFT solver may not exceed the number of Hubbard atoms, which will usually restrict us to a single task. Gains can be made using OMP threading - for the solver, at least.

Monitoring convergence
----------------------
For any given run, it is important to monitor a number of different parameters. Firstly, the goodness of the AIM mapping is quantified by :math:`d`; this ought to be small (on the order of $10^{-7}$ at least as a rough guide).

Several scripts are available to help monitor convergence. These can be found in...

For a more detailed view, it is possible to visually assess the AIM fitting procedure. This can be done by entering the
``atom<#>/dir_green_outputX_X_X_iterX`` directory and running ``check_convergence.out``. This will produce a series of files ``compare_imp_<site index>_<spin index>``: these contain the real and imaginary parts of the impurity Green function :math`G_{loc}(\omega)`\ . There is an equivalent set ``compare_lat`` that contain the local projected Green function :math`$G_{imp}(\omega)`\ . If the solver well parameterised and working, these two Greens functions will match.

These plots can help inform a good choice othe parameters `nw`, `fit_weight...`

Common sources of convergence failure
-------------------------------------
Stagnating distance parameter: provide the bath with additional degrees of freedom, for example by increasing ``sites_ed`` or ``ncpt_approx``.
Disagreement at low freqeuncies: increase the fit weight parameter to prioritise low frequencies. Capturing low-frequency behaviour accurately is crucial.

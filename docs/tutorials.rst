Tutorial
========

In this tutorial we will run a DFT+DMFT calculation on FeCO\ :sub:`3`\ .

.. raw:: html

   <div style="align: left; text-align:center;">
      <img src="_static/tutorial_files/feco3.png" width="40%">
      <div class="caption">
      The FeCO<sub>3</sub> molecule 
      </div>
   </div>


The files for this tutorial can be found in the TOSCAM git repository under ``tutorials/``

DFT geometry optimisation
-------------------------

DMFT
----

Removal of CPT approximation
----------------------------

Density of States
-----------------

Optical spectra
---------------
Optical spectra can be calculated as an additional post-processing step. To understand what's going involved in this calculation, perform the following steps:

1. In the ``<seed>.dat`` file, alter the following settings
::

   dmft_optics : T
   dmft_optics_i1 : 1
   dmft_optics_i2 : 1

2. Run a standard ONETEP calculation with ``calulation: properties``. This will produce, among other things, a ``store_nabla1`` file
3. run ``onetep.dmft.compute.dos.out``
4. run ``onetep.dmft.compute.optics.out``. The optical conductivity :math:`\sigma_{xx}` is then contained in ``optical_conductivity_pm/spin1/spin2_vertex_1``. To calculate the optical conductivity along other axes, change ``dmft_optics_i1`` and ``i2``; 1, 2, and 3 correspond to the :math:`x`, :math:`y`, and :math:`z` axes respectively.

The isotropic optical conductivity :math:`\frac{1}{3}\left(\sigma_{xx} + \sigma_{yy} + \sigma_{zz}\right)` can be straightforwardly calculated in a single step using the shell script `isotropic_anisotropic_optical_tensor.out`
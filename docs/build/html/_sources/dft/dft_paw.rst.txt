.. _dft_paw:

==================================================
Projector Augmented Wave Density Functional Theory
==================================================

Let us motivate our discussion of the specifics of projector augmented wave density functional theory (PAW-DFT):cite:`blochl1994_paw,kresse1999_paw`, by revisiting the implication of crystalline arrangements of atoms on the electron density.  Suppose we have the real space lattice vectors: :math:`\mathbf{a}_1`, :math:`\mathbf{a}_2`, and :math:`\mathbf{a}_3`, which produces the Bravais lattice :math:`\mathbf{R}_{\mathbf{n}}`, where :math:`\mathbf{n}=\left[n_1,n_2,n_3\right]`.

.. math::

   \mathbf{R}_{\mathbf{n}} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + n_3 \mathbf{a}_3

Due to periodicity, all functions dependent upon the Bravais lattice must have a periodic representation.  This includes the electronic density in an atomic crystal :math:`\rho(\mathbf{r})` which can be written as a periodic function

.. math::

   \rho(\mathbf{r}) = \rho(\mathbf{R}_{\mathbf{n}} + \mathbf{r})

And it is useful to represent it as a Fourier series expansion, an expansion in sines and cosines.

.. math::

   \rho(\mathbf{r}) = \sum_m \rho_m e^{i \mathbf{G}_m \cdot \mathbf{r}} + e^{i \mathbf{G}_m \cdot \mathbf{R}_n}

Since the function :math:`\rho(\mathbf{r})` is periodic, then for any choice of :math:`\mathbf{n},\mathbf{k} \in \mathbb{R}`


    
Then the reciprocal lattice vectors are defined


Bloch's Theorem
---------------
The electron wavefunctions of a crystal have a basis consisting entirely of Bloch wave energy eigenstates.

The energy eigenstates for an electron in a crystal can be written as Bloch waves.

A wavefunction :math:`\Psi` is a Bloch wave if it has the form:

.. math::

   \Psi(\mathbf{r}) = e^{ik \cdot r} u(\mathbf{r})

Energy Cutoff
-------------

For a DFT calculation with plane waves, the electronic wavefunction is represented as the infinite summation of plane waves, which mus be truncated to a finite series.

With more plane waves, there is more accuracy, but also at higher computational cost.  

Plane waves with less kinetic energy

.. math::

   \frac{\hbar^2}{2m}\lvert \bm{k} + \bm{G} \rvert^2

have a higher contribution to the sum, so the plane waves with lower energy have the highest contribution.

the deermination of the energy cutoff.  
Detemining of the energy cutoff for the plane wave basis set expansion has a large effect on the cost of calculation as well as the accuracy of calculation.  If :math:`E_{cut}` is the energy cutoff, then plane waves with kinetic energy less than :math:`E_{cut}` are excluded from the basis set.

.. math::

   \lvert \mathbf{G}+\mathbf{k} \rvert < G_{cut}

where

.. math::

   E_{cut}=\frac{\hbar}{2m}G^2_{cut}

Additional Reading
==================

* `Projector Augmented-Wave Method (PAW) <https://github.com/certik/sphinx-jax/blob/master/src/quantum/paw.rst>`_

* `The Projector Augmented-wave method <https://arxiv.org/pdf/0910.1921v2.pdf>_ by Carsten Rostgaard
References
----------

.. bibliography:: dft_paw.bib

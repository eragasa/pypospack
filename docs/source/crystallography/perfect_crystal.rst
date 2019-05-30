.. _perfect_crystal:
===========================================
The Perfect Crystal and the Simulation Cell
===========================================

The Perfect Crystal
===================

Inorganic materials are typically crystalline, meaning that they are periodic at the atomic scale.  
A crystal consists of atoms arranged in a pattern that repeats periodically in three dimensions.  In this defiinition, the pattern can consist of a single atom, a group of atoms, a molecule or group of molecules.  As a solid forms, atoms and/or molecules asume fixed orientations and positions with respect to each other.  As a necessary consequence of particle growth, atoms and molecules position themselves as to minimize the forces acting upon it.  Each molecule entering the solid phase is influenced in almost exactly the same way as the proceeding molecule, and the solid particle consists of three-dimensional ordered arrray of molecules; that is a crystal.

A crystal is a periodic array of atoms and consists of a set of lattice vectors which all of Euclidean space and an atomic basis, which when combined form a unit cell.  For computational materials, it is convenience to represent crystals as a material with infinite extent, but with a finite representation for computational tractibility.  The coordinate system in which to represent an infinite solid.

See books by 

:cite:`sands1969_book_crystals` Donald Sands. `Introduction to Crystallography <https://www.amazon.com/Introduction-Crystallography-Dover-Books-Chemistry-ebook/dp/B008TVLYUC/ref=sr_1_2?ie=UTF8&qid=1502437192&sr=8-2&keywords=crystallography>`_.  This book is currently published by Dover publications and is a rather cheap book on the subject.

:cite:`sands2002_book_crystals_2`. Donald Sands. `Vectors and Tensors in Crystallography <https://www.amazon.com/Vectors-Tensors-Crystallography-Donald-Sands/dp/0201071479/ref=sr_1_1?ie=UTF8&qid=1502437679&sr=8-1&keywords=vectors+and+tensors+in+crystallography>`_.  This is a much more mathematical treatment of crystallography. 

:cite:`hammond2001_book_crystal` Christopher Hammond. `Basics of Crystallography <https://www.amazon.com/Basics-Crystallography-Diffraction-Fourth-International/dp/0198738684/ref=sr_1_5?ie=UTF8&qid=1502436584&sr=8-5&keywords=crystallography>`_


:cite:`giacovazzo2002_book_crystal` Giacovazzo and Monaco. `Fundamentals of Crystallography <https://www.amazon.com/Fundamentals-Crystallography-C-Giacovazzo/dp/0198509588/ref=sr_1_3?s=books&ie=UTF8&qid=1502436703&sr=1-3&keywords=giacovazzo>`_

A Graphical Explanation of a Lattice
------------------------------------

We will replace the traditional terminology of the unit cell with a mathematical description which is more suitable for computational materials.  The notation for computational materials isn't standardized, but the terminology selected here is consistent with computational materials research, and the exposition here is one for practical notation and amenable for the application of applied mathematics and understanding computational underpinings.

In conventional notation, the lattice vectors are often referred to as :math:`\mathbf{a}`, :math:`\mathbf{b}`, and :math:`\mathbf{c}`, and are often conveniently described with the length of the lattice vectors, :math:`a`, :math:`b`, and :math:`c`, to represent the length of the vectors :math:`\mathbf{a}`, :math:`\mathbf{b}`, and :math:`\mathbf{c}`, respectively.  The angles :math:`\alpha`, :math:`\beta`, and :math:`\gamma` are used to describe the angles between the lattice vectors.

A Mathematical Explanation of a Lattice
---------------------------------------

Computationally, these are an inelegant representation of a crystal lattice and use the following notation instead.  These lattice vectors are denoted :math:`\mathbf{a}_1`, :math:`\mathbf{a}_2`, and :math:`\mathbf{a}_3`.  These vectors are often denoted as the row vectors :math:`\mathbf{a}_{i} = [a_{i1}, a_{i2}, a_{i3}]`, but are often more useful in matrix operations as column vectors.  Whether or not the lattice vectors should be used as a column or row vector, should be obvious by the application.

In computational materials, these lattice vectors are often represented as the matrix :math:`\mathbf{H}`, where

.. math::

   \mathbf{H} 
   = 
   \begin{bmatrix}
       a_{11} & a_{12} & a_{13} \\
       a_{21} & a_{22} & a_{23} \\
       a_{31} & a_{32} & a_{33}
   \end{bmatrix}
   =
   \begin{bmatrix} 
       \mathbf{a}_1 \\
       \mathbf{a}_2 \\ 
       \mathbf{a}_3 
   \end{bmatrix}

Cartesian Coordinate System and Direct Lattice Coordinates
----------------------------------------------------------

There are two types of representation of the basis of atoms in computational material system.  One represents the basis of atoms in the unit cell using the standard unit vectors in Eucliean space, which is commonly known as the Cartesian coordinate sytem.  The second in the representation of the lattice vectors of the atom.

Given the representation of :math:`\mathbf{H}` as a matrix containing the lattice vectors, then the transpose of that matrix is :math:`\mathbf{H}^T`.  The transformation of the direct coordinates to cartesian coordinates is 

.. math::
   
   \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \mathbf{H}^T \begin{bmatrix} u \\ v \\ w \end{bmatrix}

where :math:`x`, :math:`y`, and :math:`z` are elements of a point in a Cartesian coordinate system and :math:`u`, :math:`v`, and :math:`w` are elements of the same points defined by the edge vectors of the lattice.

.. math::
   \begin{bmatrix} 
        x\\ y \\ z 
   \end{bmatrix} 
   =
   \begin{bmatrix}
       a_{11} & a_{21} & a_{31} \\
       a_{12} & a_{22} & a_{32} \\
       a_{13} & a_{23} & a_{33}
   \end{bmatrix}
   \begin{bmatrix} 
       u \\ v \\ w 
   \end{bmatrix}

The transformation from the Cartesian coordinate system to the direct lattice is

.. math::

   \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \mathbf{H}^T \begin{bmatrix} u \\ v \\ w \end{bmatrix}

   \mathbf{H}^T)^{-1} \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{matrix} u \\ v \\ w \end{matrix}

The Simulation Cell
===================

The simulation cell consists of the lattice vectors and a basis of atoms in that lattice vector.  Atoms of different species must be identified.  From a notational perspective, differentiation of the chemical species uses small Greek subscripts to identify different species as a general terminology (e.g. :math:`\mathbf{r}_{\alpha}` or the chemical symbols for specificity (e.g. :math:`\mathbf{r}_{Cu}`).

The atomic basis consists of the atoms which are contained within the volume bounded by the lattice vectors.  From a programatic implementation, these basis atoms are implemented as an Atom class which as the properties: symbol, position, index, and variety of specific properties (such as magnetic moment).

Lattice Vectors and Periodic Boundary Conditions
------------------------------------------------

In order to motivate the idea of an infinite solid mathematically, let us start with a function which is periodic in one dimension.  A function :math:`f` is periodic with period :math:`P` if

.. math::

   f(x+P) = f(x)

where the periodic function can be defined as a function whose graph exhibits translational symmetry.  The most ubiqitious of these function are the trigonmetric series which are integral to calculating and modelling solid matter.  The sine function is periodic

.. math::

   \sin(x+2\pi) = \sin(x)

This concept of periodicity is also referred to as translational invariance.  The set of points due to translational symmetry from the infinite discrete set

.. math::

   \{ x + n P \}

In three-dimensions, translational symmetry are in the direction of lattice vectors.  If periodic boundary conditions are applied in three directions, the three lattice vectors must span all of Euclidean 3D space.  

We then motivate the discussion through molecular dynamics and the strength in molecular dynamics in solving dynamical systems by sampling of the Boltzmann distribution through a variety of thermodynamic ensembles by looking at the Carr-Parrinello method.  Due to the computational cost of a DFT calculation between the updates in the Carr-Parrinello method, the need for empirical interatomic potentials becomes motivated.  We discuss different classes of interatomic potentials.

Additional References
=====================

References
=========
.. bibliography:: perfect_crystal.bib

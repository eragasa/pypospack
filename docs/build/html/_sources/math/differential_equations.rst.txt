===================================================
An Computational Apprach to Differnential Equations
===================================================


Existance and Uniqueness of Solutions
=====================================

The existence and uniqueness theorems in differential equations helps us determine if a solution to differential equation and if that solution is unique.

To motivate this example let us consider existance and uniqueness of solutions for simple algebraic systems.  Consider :math:`3x+4=0` has a solution for :math:`x` exists, :math:`x=4/3`, which is also a unique equation.  However, the quadtratic equation :math:`x^2=4` has the solution :math:`x=\{-2,2\}`.  In this situation, the quadratic equatihas solutions that exist, but those solutions are not unique.

Existance and Uniqueness Theorem.  Given the initial boundary problem (IVP)

.. math::
  
   \frac{dy}{dt} = f(t,y), y(t_0)=y_0

(1) Exitance.  If :math:`f(t,y)` is continuous in the region surrounding the initial value, :math:`y(t_0)=y-0`, then we can define a rectangular region :math:`R`, 

.. math::

   R = \bigl \{(t,y)|a<t<b,c<y<d \bigr \}

that contains the point :math:`(t_0,y_0)`, then the IVP has a unique solution, y(t)

(2) Uniqueness.  If :math:`\partial f / \partial y` is continuous in the region :math:`R` then the IVP has a unique solution, y(t)

=================
Fourier Transform
=================

James O'Brien "`Introduction to Fourier Transforms for Physicists. <http://www.phys.uconn.edu/~obrien/index_files/fourier.pdf>`_"

Let us denote, the Fourier transform of the function :math:`f` is denoted by the hat symbol: :math:`\hat{f}`.

.. math::

   \hat{f}(\xi) = \int_{-\infty}^{\infty} f(x) exp(-2 \pi i x \xi) dx

for any real number :math:`\xi`

Under suitable conditions, :math:`f` can be determined from :math:`\hat{f}` through the inverse transform

.. math::

   f(x) = \int_{\infty}^{\infty} \hat{f}(\xi) exp(-2 \pi i x \xi) d \xi


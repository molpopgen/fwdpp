.. _rng:

Random number generators
+++++++++++++++++++++++++++++++++++++

The random number generator type is :type:`fwdpp::GSLrng_mt`, which is a mersenne twister.
The class constructor accepts a single, 32-bit unsigned integer.
This class is effectively a ``unique_ptr`` with a custom deleter.
As such, it defines a move-only type.

:type:`fwdpp::GSLrng_mt` is defined as a template ``typedef`` of :class:`fwdpp::GSLrng_t`.
The function :func:`fwdpp::GSLrng_t::get` returns the underlying `const gsl_rng *`.

These types are defined in :ref:`file_fwdpp_GSLrng_t.hpp`.

To use this type:

.. literalinclude:: ../../examples/rng_example.cc

.. note::

   Most or all of ``fwdpp`` is compatible with using a "bare" ``gsl_rng *``.
   However, you give up exception safety by doing so.

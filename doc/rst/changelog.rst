.. _changelog:

Change log
=============================================

0.9.0
********************************************

This is a pretty big release.
In fact, it is too big to list all of the changes here.
They are collected on GitHub under the ``0.9.0`` milestone.

Importantly, this release deprecates several features.
Where possible, a ``C++14``-style ``[[deprecated]]`` marker is used.

We plan to start releasing more often, so that individual releases aren't so big.

In future releases, we plan to improve the organization of the library headers.
We will also clean up the test suite to reflect the new layout.
We will keep the existing headers and use the preprocessor to emit a warning at compile time.

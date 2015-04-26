# Tutorial 5: the "sugar" layer on top of fwdpp

Through __fwdpp__ 0.2.9, a lot of the basic setup of a simulation could get tedious due to having to define a large number of container types.  In version 0.3.0 of the library, I introduced the subdirectory fwdpp/sugar, which is a layer of "syntatic sugar" designed to make fwdpp both easier to use and easier to wrap for use in other programming environments.

These functions are only relevant to _indivdual-based_ simulations.  The gamete-based portion of __fwdpp__ is unlikely to be supported (which is fine, as it is slower...).

## Overview

The sugar layer is comprised of several different headers, all located in _fwdpp/sugar_.  These header files are intended to be used directly for implementing simulations.

Further subdirectories contain the details of the namespace KTfwd::sugar, which most library users are unlikely to need to deal with directly.

Importantly, _the use of these headers is optional_, and may not work for every type of simulation.  Much of this is based upon what I've found to work for my own research using __fwdpp__.


## Smart wrapper to GSL random number generators

## Simplifying the declaration of a population

## Simplifying serializing of simulated data

## Built-in "standard" mutation types

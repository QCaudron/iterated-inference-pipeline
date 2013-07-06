Iterated inference pipeline
===========================

_Pipable_ plug-and-play inference methods for time series analysis with state space models.

    cat theta.json | ./simplex -M 10000 -P | ./ksimplex -M 10000 -P > mle.json
    cat mle.json | ./kmcmc -M 100000 --full -P | ./pmcmc -J 1000 -M 500000 --full


This README provided information for developers or users wanting to
build pipable inference methods on their local machine. If you want
to **use** those methods, see [plom.io](http://plom.io/cli)

[plom.io](http://plom.io/cli) provides a compilation service that
generate standalone binaries that you can run on your machine.

Overview
========

All the methods are implemented in plain C.  The C code contain
generic part (working with any models) and model specific part.  The
specific parts are templated using Python.

Installation
============

##Dependencies

C:
- gsl: http://www.gnu.org/software/gsl/
- zmq: http://www.zeromq.org/
- jansson: http://www.digip.org/jansson/
- openMP: http://openmp.org/ (openMP is optional as it is not supported by Clang. All the feature provided by openMP have a counterpart provided by zmq and pthread.)

Python:
- Python 2.7.x: www.python.org/
- Django: https://www.djangoproject.com/
- SymPy: http://sympy.org/

##Building the C libraries
in model_builder/C:

    make
    make install

##Creating and installing the python package

At the root of the repo run:

    python setup.py sdist

in the package directory:

    python setup.py install


Usage
=====

##Generating the model-specific code:

From the command line, run:

    pmbuilder context.json process.json link.json -o path_model_coded_in_C

In your script you can use:

    import json
    from plom.Builder import PlomModelBuilder

    c = json.load(open('example/noise/context.json'))
    p = json.load(open('example/noise/process.json'))
    l = json.load(open('example/noise/link.json'))

    model = PlomModelBuilder('path_model_coded_in_C', c, p, l)

    model.prepare()
    model.write_settings()
    model.render()

##Building the inference methods

in path_model_coded_in_C/C/templates:

    make
    make install
    
All the inference methods binaries are now available in
path_model_coded_in_C


Contributing
============

<!--
First generate the documentation for the C code:
A general introduction:
In model_builder/doc/ run: ```docco smc.c``` and open docs/smc.html with
a web browser. 
In the same vain, ```docco kalman.c``` will detail our implementation
of the extended Kalman Filter.
-->

On model_builder/doc/ run: ```doxygen Doxyfile``` and open
model_builder/doc/html/index.html with a web browser.

After having learned the basic structures involved in ```core```, we
recommend to use call and caller graphs as an entry point to the
sources.

Tests
=====

In model_builder:

    python -m unittest discover


License
=======

GPL version 3 or any later version.


Acknowledgements
================

We want to thank:

- Professor
  [Bryan Grenfell](http://www.princeton.edu/eeb/people/display_person.xml?netid=grenfell).

- Professor
  [Bernard Cazelles](http://www.biologie.ens.fr/~cazelles/bernard/Welcome.html).

- Princeton University and the "simforence" project for letting us
  release the sources under the GPLv3 license.
  

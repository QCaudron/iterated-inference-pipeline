plom-sfi
========

**Plug-and-play** inference methods in plain C for http://www.plom.io/.

usage
=====

This README provided information for developers
or users wanting to build plom-sfi on their local machine.

If you want to use the method provided by plom-sfi see
http://www.plom.io/cli

You do not need to install plom-sfi or its dependencies to use
it. plom.io provide a compilation service that generate standalone
binaries that you can run on your machine so go to
http://www.plom.io/cli

introduction
============

The C code contain generic part and model specific part. The C code of
the specific parts are templated using the plom python package.


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

    ./install.sh

(see source for details)


##Usage

###Generating the model-specific code:

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

###Building the inference methods

in path_model_coded_in_C/C/templates:

    make
    make install
    
All the inference methods binaries are now available in path_model_coded_in_C

See http://www.plom.io/cli for documentation on how to
use the generated binaries.


##Contributing to the C library

First generate the documentation for the C code:
A general introduction:
In model_builder/doc/ run: ```docco smc.c``` and open docs/smc.html with
a web browser. 
In the same vain, ```docco kalman.c``` will detail our implementation
of the extended Kalman Filter.

in model_builder/doc/ run: ```doxygen Doxyfile``` and open model_builder/doc/html/index.html with
a web browser.

After having learned the basic structures involved in ```core```, we
recommend to use call and caller graphs as an entry point to the
sources.

##Tests

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

- Princeton University for letting us release the code of plom-sfi
  under the GPLv3 license.

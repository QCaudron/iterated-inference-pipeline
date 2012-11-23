plom-sfi
========

Simulation Forecasting and Inference methods for PLoM

NOTE: This README is dedicated to Plom developers.  If you are
interested in using Plom please go to http://www.plom.io/
where you will find appropriate documentation.

- plom command line tools. This code lives in script/
- Python code to generate model in plain C. This code lives in model_builder/
- C code to run models. This code lives in model_builder/C_lib/src/C/ The C code contain generic part and model specific part that are rendered using model_builder

##Dependencies

C:
- gsl: http://www.gnu.org/software/gsl/
- zmq: http://www.zeromq.org/
- jansson: http://www.digip.org/jansson/
- openMP: http://openmp.org/

Python:
- Python 2.7.x: www.python.org/
- Django: https://www.djangoproject.com/
- SymPy: http://sympy.org/
- NumPy: http://numpy.scipy.org/


##Creating and installing plom python package (containing the C code as package data)

At the root of the repo run:

    ./install.sh

(see source for details)


##Usage

See http://www.plom.io/doc/modeler/intro


##Contributing to the C library

First generate the documentation for the C code:
in model_builder/C_lib/doc/ run: ```doxygen Doxyfile``` and open model_builder/C_lib/doc/html/index.html with
a web browser.

After having learned the basic structures involved in ```core```, we
recommend to use call and caller graphs as an entry point to the
sources.


License
=======

GPL version 3 or any later version.

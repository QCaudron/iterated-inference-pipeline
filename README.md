plom-sfi
========

Simulation Forecasting and Inference methods for PLoM

NOTE: This README is dedicated to plom-sfi developers.  If you are
interested in using PLoM simulation, forecasting and inference methods
please go to http://www.plom.io/ where you will find appropriate
documentation.

- pmbuilder command line tools. The code lives in script/
- Python code to generate model in plain C. The code lives in model_builder/
- C code to perform simulation, forecasting and inference. This code lives in model_builder/C_lib/src/C/ The C code contain generic part and model specific part that are rendered using model_builder

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


##Creating and installing the python package (containing the C code as package data)

At the root of the repo run:

    ./install.sh

(see source for details)


##Usage

See http://www.plom.io/doc/modeler/intro


    import json
    from plom.Builder import PlomModelBuilder

    c = json.load(open('context.json'))
    p = json.load(open('process.json'))
    l = json.load(open('link.json'))

    model = PlomModelBuilder('path_model_coded_in_C', c, p, l)

    model.prepare()
    model.write_settings()
    model.code()
    model.compile()



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

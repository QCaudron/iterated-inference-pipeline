#!/usr/bin/env python

from distutils.core import setup

setup(name='plom',
      version='0.2.0',
      description='PLoM model builder for epidemiology',
      author='Sebastien Ballesteros',
      author_email='sebastien@plom.io',
      url='http://www.plom.io',
      packages=['plom'],
      package_dir={'plom': 'model_builder'},
      scripts=['scripts/pmbuilder'],
      package_data={'plom': ['C/lib/*',
                             'C/core/*',
                             'C/templates/*',
                             'C/kalman/*',
                             'C/mif/*',
                             'C/pmcmc/*',
                             'C/simplex/*',
                             'C/simulation/*',
                             'C/smc/*',
                             'C/worker/*']}
)

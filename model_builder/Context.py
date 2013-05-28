##########################################################################
#    This file is part of plom.
#
#    plom is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    plom is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with plom.  If not, see
#    <http://www.gnu.org/licenses/>.
#########################################################################

import csv
import copy
import json
import sys
import os
#from dateutil import rrule
import datetime

class PlomContextError(Exception):
    pass


class Context:
    """handles a Data Block"""

    def __init__(self, context, **kwargs):
        """Create a Context object
        """

        self.ts_id = [x['id'] for x in context['time_series']] ##NOTE: we can't guarantee an order on ts_id until it is linked with a model
        self.map_ts_cac = {x['id']: x['population_id'][:] for x in context['time_series']}

        self.cac_id = [x['id'] for x in context['population']]


        _headerline = [x.split('__') for x in self.cac_id]
        _repeated_cities = [ x[0] for x in _headerline ]
        _repeated_ages = [ x[1] for x in _headerline ]
        self.cities_id = [ _repeated_cities[x] for x in sorted([ _repeated_cities.index(y) for y in set(_repeated_cities) ]) ]
        self.ages_id = [ _repeated_ages[x] for x in sorted([ _repeated_ages.index(y) for y in set(_repeated_ages) ]) ]
        self.N_C = len(set(_repeated_cities))
        self.N_AC = len(set([x[1] for x in _headerline]))
        self.N_CAC = self.N_C*self.N_AC


        _headerline = [x.split('__') for x in self.ts_id]
        self.N_TS = len(_headerline)
        self._repeated_name_ts = [x[0] for x in _headerline]
        self._repeated_name_stream = [x[1] for x in _headerline]
        self._repeated_obs_type = [x[2] for x in _headerline]

        self.N_TS_INC = len( [x for i, x in enumerate(self._repeated_name_ts) if self._repeated_obs_type[i] == 'inc'] )
        self.N_TS_INC_UNIQUE = len(set( [x for i, x in enumerate(self._repeated_name_ts) if self._repeated_obs_type[i] == 'inc'] ))


        ##############################################
        ##data
        ##############################################
        self.par_fixed_values = {}

        ##First we ensure that mandatory properties are represented in self
        self.data = []; self.dates = []; self.prop = []

        for d in (context.get('data', []) + context.get('metadata', [])):
            if d['id'] == 'data':
                mydata = self.handle_context_data(d['source'])
                self.data = mydata['values']
                self.dates = mydata['dates']
            elif d['id'] == 'prop':
                self.prop = self.handle_context_data(d['source'])['values']
            else:
                self.par_fixed_values[d['id']] = self.handle_context_data(d['source'])['values']

        self.N_DATA = len(self.data)

        self.N_DATA_PAR_FIXED = len(self.par_fixed_values[self.par_fixed_values.keys()[0]]) if self.par_fixed_values else 1
        if not self.prop:
            print('\033[93m' + "WARNING: " + '\033[0m' + 'property "prop" (proportion under surveillance) is not present in context.data, default to 1.0')
            self.prop = [[1.0]*self.N_TS]*self.N_DATA_PAR_FIXED

        ##TO DO:
        ##self.school_terms = copy.deepcopy(school_terms)


        self.frequency = context.get('frequency', None) ##determined by plom(1) (recommended way)
        if not self.frequency: ##quick and dirty
            delta = datetime.datetime.strptime(self.dates[1], "%Y-%m-%d").date() - datetime.datetime.strptime(self.dates[0], "%Y-%m-%d").date()
            self.frequency = {1: 'D', 7: 'W', 30: 'M', 365: 'Y'}[delta.days]


    def handle_context_data(self, source):
        """
        Takes into account that source can be a path to a csv or a native array and check that the header match the context header.
        """

        #if source is a string : read from csv and replace source by the read array
        if isinstance(source, str) or isinstance(source, unicode):
            try :
                f = open(source, 'r')
            except IOError:
                raise PlomContextError('\033[91m' + 'FAILURE! ' + '\033[0m' + source + ' from ' + os.getcwd() + ' could not be found')
            else:
                reader = csv.reader(f, delimiter=',', quotechar='"')
                header = [next(reader)] #skip header
                source = header + [[row[0]] + map(lambda x: float(x) if x!='null' else None, row[1:]) for row in reader if row!=[]]
                f.close()

        #source is now an array (it was or it has been transformed into one with csv)
        if source:
            if source[0][0] != 'date':
                raise PlomContextError("\033[91mFAIL:\033[0m data header doesn't match context description {0} != date".format(source[0][0]))
            if source[0][1:] != self.cac_id and source[0][1:] !=self.ts_id:
                raise PlomContextError("\033[91mFAIL:\033[0m data header doesn't match context description")

            data = {'header': source[0],
                    'values': [x[1:] for x in source[1:]],
                    'dates': [x[0] for x in source[1:]]}

        else:
            data = {'header':[], 'values':[], 'dates':[]}

        return data


if __name__ == '__main__':

    c = json.load(open(os.path.join('example', 'noise','context.json')))
    ##fix path (this is normally done by pmbuilder(1))
    for x in c['data'] + c['metadata']:
        x['source'] = os.path.join('example', 'noise', x['source'])

    context = Context(c)

##    print(context.prop)
##    print(context.cac_id)
##    print(context.par_fixed_values)
##    print(context.data)

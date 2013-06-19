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

import copy
import json
import sys
import os
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
        self.data = []; self.dates = []

        if 'source' in context['data']:            
            self.dates = [x[0] for x in context['data']['source'][self.ts_id[0]]['value']]

            for i in range(len(self.dates)):
                self.data.append([context['data']['source'][ts]['value'][i][1] for ts in self.ts_id])
            
        for d in context.get('metadata', []):        
            self.par_fixed_values[d['id']] = copy.deepcopy(d)

        self.N_DATA = len(self.data)
        
        self.date_0 = datetime.datetime.strptime(context['data']['t0'], "%Y-%m-%d").date()


        ##TO DO:
        ##self.school_terms = copy.deepcopy(school_terms)




if __name__ == '__main__':

    c = json.load(open(os.path.join('example', 'noise','context.json')))

    context = Context(c)

    print(context.cac_id)
    print(context.data)
    print(context.par_fixed_values)

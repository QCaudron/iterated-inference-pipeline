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

import os
import os.path
import json
import datetime
import copy

def get_obs2ts(self):
    """get obs2ts: list with for every obs_var:
    -n_ts_unique
    for every n_ts_unique:
    -n_stream
    -n_cac
    -cac (list of tuples containing index of cities and age classes aggregated in the time serie)

    NOTE: this function assume that data have been sorted by self.obs_var and stream
    """

    obs2ts = []

    for o in self.obs_var:

        ind_ts_o = [self.ts_id.index(x) for x in self.ts_id if self.map_ts_obs[x] == o]
        ts_o = [self._repeated_name_ts[x] for x in ind_ts_o]

        ts_o_unique = [] ##we can't use set(ts_o) as we need to preserve order
        for e in ts_o:
            if e not in ts_o_unique:
                ts_o_unique.append(e)

        n_ts_o_unique = len(ts_o_unique)

        info_o = {}
        info_o['n_ts_unique'] = n_ts_o_unique
        info_o['n_stream'] = []
        info_o['n_cac'] = []
        info_o['cac'] = []

        for ts_o_u in ts_o_unique:
            ind_ts_o_u = [ ind_ts_o[i] for i, x in enumerate(ts_o) if x == ts_o_u ]
            stream_ts_o_u = [ self._repeated_name_stream[x] for x in ind_ts_o_u ]
            n_stream_ts_o_u = len(set(stream_ts_o_u))

            cac_ts_o_u = self.map_ts_cac[ self.ts_id[ ind_ts_o_u[0] ] ] ##all ts from ind_ts_o_u have the same cac
            n_cac_ts_o_u = len(cac_ts_o_u)

            info_o['n_stream'].append(n_stream_ts_o_u)
            info_o['n_cac'].append(n_cac_ts_o_u)
            info_o['cac'].append( [(self.cities_id.index(x.split('__')[0]), self.ages_id.index(x.split('__')[1]))  for x in cac_ts_o_u]  )

        obs2ts.append(info_o)

    ##add offset: first index of ts for an observed variable
    offset=0
    for x in obs2ts:
        x['offset'] = offset
        offset += sum(x['n_stream'])

    return obs2ts



def make_settings_json(self):

    settings = {}

    settings['POP_SIZE_EQ_SUM_SV'] = (self.remainder == None)

    #######data
    settings['data'] = {}

    ##data/obs2ts
    settings['data']['obs2ts'] = get_obs2ts(self)

    ##data/data (be sure to have sorted the context before this part)
    settings['data']['data'] = self.data
    settings['data']['par_fixed_values'] = copy.deepcopy(self.par_fixed_values)

    #convert dates into days since first data points
    for k, v in settings['data']['par_fixed_values'].iteritems():
        times = []
        for i, d in enumerate(v['dates']):
            delta =  datetime.datetime.strptime(d, "%Y-%m-%d").date() - self.date_0
            times.append(delta.days)

        #change values format so that it's easy to generate GSL interpolators'
        values = []        
        for i in range(len(v['values'][0])):                    
            obj = {'id': v['header'][i+1], 'x': times, 'y': [v['values'][x][i] for x in range(len(v['values']))], 'size': len(v['values']) }
            if 'unit' in v:
                obj['unit'] = v['unit']
            if 'type' in v:
                obj['type'] = v['type']
            values.append(obj)
        
        settings['data']['par_fixed_values'][k] = values
        

    settings['data']['dates'] = self.dates
    settings['data']['times'] = [0] + [(datetime.datetime.strptime(x, "%Y-%m-%d").date() - self.date_0).days for x in self.dates]

    ##TODO
    ##settings['data']['school_terms'] = self.school_terms

    ##drift
    all_order = self.par_sv + self.par_proc + self.par_obs
    settings['drift'] = {'ind_par_Xdrift_applied': [ all_order.index(x) for x in self.drift_par_proc + self.drift_par_obs ],
                         'ind_volatility_Xdrift': [ all_order.index(x) for x in self.vol_par_proc + self.vol_par_obs ]}

    settings['ind_noise_sd'] = list(set([all_order.index(x['sd']) for x in self.white_noise])) ##set as different noise can have the same intensity
    
    settings['remainder'] = self.remainder
    settings['date_0'] = str(self.date_0)

    #######cst settings
    settings['cst'] = {'N_C': self.N_C,
                       'N_AC': self.N_AC,
                       'N_PAR_PROC': len(self.par_proc),
                       'N_PAR_OBS': len(self.par_obs),
                       'N_PAR_SV':  len(self.par_sv),
                       'N_PAR_FIXED':  len(self.par_fixed),
                       'N_TS': self.N_TS,
                       'N_TS_INC': self.N_TS_INC,
                       'N_TS_INC_UNIQUE': self.N_TS_INC_UNIQUE,
                       'N_DATA': self.N_DATA,
                       'N_OBS_ALL': len(self.obs_var_def),
                       'N_OBS_INC': len([x for x in self.obs_var_def if isinstance(x[0], dict)]),
                       'N_OBS_PREV': len([x for x in self.obs_var_def if not isinstance(x[0], dict) ]),
                       'N_DRIFT':len(self.drift_par_proc) + len(self.drift_par_obs),
                       'IS_SCHOOL_TERMS':1 if any(map(lambda x: 'terms_forcing' in x, [r['rate'] for r in self.proc_model])) else 0}

    #######order settings (be sure to have sorted the context before this part)
    settings['orders'] = {}
    settings['orders']['par_sv'] = self.par_sv
    settings['orders']['par_proc'] = self.par_proc
    settings['orders']['par_fixed'] = self.par_fixed
    settings['orders']['par_obs'] = self.par_obs
    settings['orders']['ts_id'] = self.ts_id
    settings['orders']['cac_id'] = self.cac_id
    settings['orders']['drift_var'] = self.drift_var

    return json.dumps(settings)

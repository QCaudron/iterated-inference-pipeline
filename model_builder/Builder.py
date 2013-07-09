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
import tarfile
import shutil
import subprocess
import copy
import datetime
import json
from Ccoder import Ccoder
from django.conf import settings as django_settings

from Context import Context

##to do: get rid of django and use Jinja2
from django.template import Context as DjangoContext
from django.template.loader import get_template

#from dateutil import rrule
#from datetime import datetime
##freq2rrule = {'D': rrule.DAILY, 'W': rrule.WEEKLY, 'M': rrule.MONTHLY, 'Y': rrule.YEARLY}
##self.dates = map(lambda x: x.strftime('%Y/%m/%d'), list(rrule.rrule(freq2rrule[self.frequency], count=self.N_DATA, dtstart=date_start) ))


def prepare_model(path_rendered, path_templates, replace=True):
    """
    copy templates to path_rendered
    """

    ##this function is called only when a new user has created or edited a model whose name is unique (primary key) so it is the only one able to recreate a model...
    if replace:
        if os.path.exists(path_rendered):
            shutil.rmtree(path_rendered)

    #copy templates to uploads/rendered/user_name/model_id
    if not os.path.exists(path_rendered):
        shutil.copytree(path_templates, os.path.join(path_rendered, 'C'))


def archive_model(path_rendered, replace=True):
    """make a tarball"""

    tar = tarfile.open(os.path.join(os.path.dirname(path_rendered), os.path.basename(path_rendered)+'.tar.gz'), "w:gz")
    tar.add(path_rendered, arcname=os.path.basename(path_rendered))
    tar.close()

    if replace:
        if os.path.exists(path_rendered):
            shutil.rmtree(path_rendered)



class PlomModelBuilder(Context, Ccoder):
    """Bind context to model and add link"""

    def __init__(self, path_rendered, context, process, link):


        Context.__init__(self, context)
        Ccoder.__init__(self, context, process, link)

        self.path_rendered = path_rendered

        ##map_ts_obs
        self.map_ts_obs = {}
        for x in link['observed']:
            for ts in x['time_series_id']:
                self.map_ts_obs[ts] = x['id']


        ##########################
        ##sort context
        ##########################

        #tbs: to be sorted
        tbs = zip(self.ts_id, self._repeated_name_stream, self._repeated_name_ts)
        #sort by data_stream
        tbs.sort(key=lambda x: x[1])        
        #sort by name_ts (in python, sorts are guaranteed to be stable)
        tbs.sort(key=lambda x: x[2])
        #sort by obs_var (in python, sorts are guaranteed to be stable)
        tbs.sort(key=lambda x: self.obs_var.index(self.map_ts_obs[ x[0] ]))


        #we need to sort ts_id, _repeated_name_ts, _repeated_name_stream, _repeated_obs_type, data and all the par_fixed involved in the obs_model
        ind_sorted = [ self.ts_id.index(x[0]) for x in tbs ]

        self.ts_id = [ self.ts_id[x] for x in ind_sorted ]
        self._repeated_name_ts = [self._repeated_name_ts[x] for x in ind_sorted]
        self._repeated_name_stream = [self._repeated_name_stream[x] for x in ind_sorted]
        self._repeated_obs_type = [self._repeated_obs_type[x] for x in ind_sorted]

        #data
        if self.data:
            self.data = [ [ y[x] for x in ind_sorted ] for y in self.data ]





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
        settings['data']['obs2ts'] = self.get_obs2ts()

        ##data/data (be sure to have sorted the context before this part)
        settings['data']['data'] = self.data
        settings['data']['par_fixed_values'] = copy.deepcopy(self.par_fixed_values)
    
        #convert dates into days since first data points
        for k, v in settings['data']['par_fixed_values'].iteritems():
            values = []
            order = self.ts_id if k in self.par_fixed_obs else self.cac_id #!!self.ts_id have been sorted
            for ts in order: 
                x = []
                y = []
                for d in v['source'][ts]['value']:
                    if d[1] != None:
                        delta =  datetime.datetime.strptime(d[0], "%Y-%m-%d").date() - self.date_0
                        x.append(delta.days)
                        y.append(d[1])


                obj = {'id': ts, 'x': x, 'y': y, 'size': len(y)}

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





    ##########################
    ##accessors
    ##########################
    def get_par_id(self):
        return {'par_sv': self.par_sv, 'par_proc': self.par_proc, 'par_obs': self.par_obs}

    def get_ts_id(self):
        "return the ordered ts_id"
        return self.ts_id

    def get_cac_id(self):
        return self.cac_id


    ##########################
    ##render model
    ##########################

    def prepare(self, path_templates=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'C'), replace=True):
        prepare_model(self.path_rendered, path_templates, replace)

    def render(self):
        """generate C code for MIF, Simplex, pMCMC, Kalman, simulation, ..."""

        if not django_settings.configured:
            django_settings.configure(TEMPLATE_DIRS = (self.path_rendered,), DEBUG = False, FILE_CHARSET = 'utf-8')

        #elif not django_settings.TEMPLATE_DIRS:
        django_settings.TEMPLATE_DIRS = (self.path_rendered,)

        is_drift = True if len(self.drift_var) > 0 else False

        ##methods whose results are use multiple times
        order = self.print_order()
        step_ode_sde = self.step_ode_sde()
        jac = self.jac(step_ode_sde['sf'])

        #core templates
        t= get_template(os.path.join(self.path_rendered, 'C', 'templates', 'core_template.c'))
    
        c = DjangoContext({'order':order,
                           'white_noise': self.white_noise,
                           'step_psr': self.step_psr(),
                           'psr_multinomial': self.psr_multinomial(),
                           'step_ode_sde': step_ode_sde,
                           'list_obs_prev': self.print_obs_prev(),
                           'obs_inc_step_psr': self.obs_inc_step_psr(),
                           'is_drift': is_drift,
                           'psr':self.print_build_psr(),
                           'proc_obs':self.print_like()})
        f = open(os.path.join(self.path_rendered, 'C', 'templates', 'core_tpl.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'C', 'templates', 'core_template.c'))

        #simulation templates
        t= get_template(os.path.join(self.path_rendered, 'C', 'templates', 'simulation_template.c'))
        c = DjangoContext({'order':order,
                           'jacobian':jac,
                           'step_ode_sde': step_ode_sde,
                           'is_drift': is_drift})
        f = open(os.path.join(self.path_rendered, 'C', 'templates', 'simulation_tpl.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'C', 'templates', 'simulation_template.c'))

        #kalman templates
        t= get_template(os.path.join(self.path_rendered, 'C', 'templates', 'kalman_template.c'))
        c = DjangoContext({'order':order,
                           'jacobian':jac,
                           'der_mean_proc_obs':self.der_mean_proc_obs(),
                           'der2_mean_proc_obs':self.der2_mean_proc_obs(),
                           'calc_Q': self.eval_Q(),
                           'is_drift': is_drift,
                           'step_ode_sde': step_ode_sde})
        f = open(os.path.join(self.path_rendered, 'C', 'templates', 'kalman_tpl.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'C', 'templates', 'kalman_template.c'))

    ##########################
    ##write
    ##########################

    def write_settings(self, settings_name='.settings'):
        with open(os.path.join(self.path_rendered, settings_name+'.json'), 'w') as f:
            f.write(self.make_settings_json())


    def archive(self, replace=True):
        archive_model(self.path_rendered, replace)


if __name__=="__main__":
    ##tutorial example...

    c = json.load(open(os.path.join('example', 'noise', 'context.json')))
    p = json.load(open(os.path.join('example', 'noise', 'process.json')))
    l = json.load(open(os.path.join('example', 'noise', 'link.json')))
            
    model = PlomModelBuilder(os.path.join(os.getenv("HOME"), 'plom_test_model'), c, p, l)

##    print model.par_fixed
##    print model.par_proc
##
##    print 'parameters: ', model.get_par_id()
##
##    print 'order of context elements: ', model.get_ts_id(), model.get_cac_id()
##
##    print model.ts_id
##    print model.cac_id
##    print model.obs_var
##    print model.map_ts_obs

    model.prepare()
    model.write_settings()
    model.render()

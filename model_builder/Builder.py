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
import json
from Ccoder import Ccoder
from plom_settings import make_settings_json
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
            f.write(make_settings_json(self))


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

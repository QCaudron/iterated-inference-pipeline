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


def cmd_exists(cmd):
    """check if a command exists"""

    return subprocess.call(["type", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0




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
        shutil.copytree(path_templates, path_rendered)


    #create settings directory (if needed)
    path_settings = os.path.join(path_rendered, 'settings')
    if not os.path.exists(path_settings):
        os.makedirs(path_settings)


def archive_model(path_rendered, replace=True):
    """make a tarball"""

    tar = tarfile.open(os.path.join(os.path.dirname(path_rendered), os.path.basename(path_rendered)+'.tar.gz'), "w:gz")
    tar.add(path_rendered, arcname=os.path.basename(path_rendered))
    tar.close()

    if replace:
        if os.path.exists(path_rendered):
            shutil.rmtree(path_rendered)



class Model(Context, Ccoder):
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

        self.partition = copy.deepcopy(link.get('partition', {}))

        ##########################
        ##Parameters
        ##########################
        #self.set_values(link)

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

        #we need to sort ts_id, _repeated_name_ts, _repeated_name_stream, _repeated_obs_type, data, prop and all the par_fixed involved in the obs_model

        #let's start easy: sort ts_id
        ind_sorted = [ self.ts_id.index(x[0]) for x in tbs ]

        self.ts_id = [ self.ts_id[x] for x in ind_sorted ]
        self._repeated_name_ts = [self._repeated_name_ts[x] for x in ind_sorted]
        self._repeated_name_stream = [self._repeated_name_stream[x] for x in ind_sorted]
        self._repeated_obs_type = [self._repeated_obs_type[x] for x in ind_sorted]

        #data and prop
        if self.prop:
            self.prop = [ [ y[x] for x in ind_sorted ] for y in self.prop ]
        if self.data:
            self.data = [ [ y[x] for x in ind_sorted ] for y in self.data ]

        #get par_fixed involved in the obs_model (par_fixed_to_be_sorted):
        par_fixed_to_be_sorted = []
        for x in self.obs_model:
            elements = self.change_user_input(x)
            for e in elements:
                if e in self.par_fixed:
                    par_fixed_to_be_sorted.append(e)

        #sort par_fixed_to_be_sorted
        if par_fixed_to_be_sorted:
            for p in par_fixed_to_be_sorted:
                self.par_fixed_values[p] = [ [ y[x] for x in ind_sorted ] for y in self.par_fixed_values[p] ]



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

    def prepare(self, path_templates=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'C_lib'), replace=True):

        prepare_model(self.path_rendered, path_templates, replace)


    def code(self):
        """generate C code for MIF, Simplex, pMCMC, Kalman, simulation, ..."""


        if not django_settings.configured:
            django_settings.configure(TEMPLATE_DIRS = (self.path_rendered,), DEBUG = False, FILE_CHARSET = 'utf-8')

        #elif not django_settings.TEMPLATE_DIRS:
        django_settings.TEMPLATE_DIRS = (self.path_rendered,)


        is_drift = True if len(self.drift_var) > 0 else False

        #core templates
        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'core', 'plom_template.h'))
        c = DjangoContext({'order':self.print_order()})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'core', 'plom.h'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'core', 'plom_template.h'))

        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'core', 'prediction_template.c'))
        c = DjangoContext({'gamma_noise': self.get_gamma_noise_terms(),
                     'print_prob': self.print_prob(),
                     'print_multinomial': self.print_multinomial(),
                     'print_update': self.print_update(),
                     'print_ode': self.print_ode(),
                     'list_obs_prev': self.print_obs_prev(),
                     'eq_obs_inc_markov': self.print_obs_inc_markov(),
                     'eq_obs_inc_ode': self.print_obs_inc_ode(),
                     'is_drift': is_drift})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'core', 'prediction.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'core', 'prediction_template.c'))

        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'core', 'build_template.c'))
        c = DjangoContext({'buildmarkov':self.print_build_markov()})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'core', 'build.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'core', 'build_template.c'))

        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'core', 'likelihood_template.c'))
        c = DjangoContext({'proc_obs':self.print_like()})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'core', 'likelihood.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'core', 'likelihood_template.c'))

        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'core', 'observation_template.c'))
        c = DjangoContext({'proc_obs':self.print_like()})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'core', 'observation.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'core', 'observation_template.c'))

        #kalman templates
        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'kalman', 'eval_ekf_update_mats_template.c'))
        c = DjangoContext({'jacobian':self.jac(),
                     'jac_proc_obs':self.jac_proc_obs,
                     'noise_Q': self.eval_Q(),
                     'stoichiometric':self.stoichiometric(),
                     'is_drift': is_drift})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'kalman', 'eval_ekf_update_mats.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'kalman', 'eval_ekf_update_mats_template.c'))

        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'kalman', 'prediction_template.c'))
        c = DjangoContext({'print_ode': self.print_ode(),
                     'eq_obs_inc_ode': self.print_obs_inc_ode(),
                     'is_drift': is_drift})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'kalman', 'prediction.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'kalman', 'prediction_template.c'))

        #simulation templates
        t= get_template(os.path.join(self.path_rendered, 'src', 'C', 'simulation', 'lyap_template.c'))
        c = DjangoContext({'jacobian':self.jac(),
                     'print_ode': self.print_ode(),
                     'is_drift': is_drift})
        f = open(os.path.join(self.path_rendered, 'src', 'C', 'simulation', 'lyap.c'),'w')
        f.write(t.render(c))
        f.close()
        os.remove(os.path.join(self.path_rendered, 'src', 'C', 'simulation', 'lyap_template.c'))


    def compile(self, web=False, simulation_only=False):
        """compiles the generated code.
        web is a flag indicating if outputs should be printed in JSON on stdout or in csv on FILE
        """

        path_src_C = os.path.join(self.path_rendered, 'src','C')

        dirs=['core', 'simulation']
        if not simulation_only:
            dirs.extend(['smc', 'simplex', 'mif', 'pmcmc', 'worker', 'kalman'])

        if web:
            cmd = 'sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 1/" plom.h'
        else:
            cmd = 'sed -ie "s/#define FLAG_JSON \([0-9]*\)/#define FLAG_JSON 0/" plom.h'

        wd = os.getcwd()

        ##set correct FLAG_JSON
        os.chdir(os.path.join(path_src_C, 'core'))

        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
        for line in p.readlines():
            print(line.rstrip())
        os.chdir(wd)

        ##compile

        cmd = 'make {0}; make install'.format('CC=gcc-4.7' if cmd_exists('gcc-4.7') else '')

        for d in dirs:
            os.chdir(os.path.join(path_src_C,  d))
            print('\033[94m  compiling : {0} (running {1})...\033[0m'.format(d, cmd))

            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
            for line in p.readlines():
                print line.rstrip()


        os.chdir(wd)



    ##########################
    ##write
    ##########################

    def write_settings(self, settings_name='settings'):
        with open(os.path.join(self.path_rendered, 'settings', settings_name+'.json'), 'w') as f:
            f.write(make_settings_json(self))


    def archive(self, replace=True):
        archive_model(self.path_rendered, replace)



if __name__=="__main__":
    ##tutorial example...

    c = json.load(open(os.path.join('examples', 'tutorial', 'context.json')))
    p = json.load(open(os.path.join('examples', 'tutorial', 'process.json')))
    l = json.load(open(os.path.join('examples', 'tutorial', 'link.json')))

    ##fix path (this is normally done by plom)
    for x in c['data']:
        x['source'] = os.path.join('examples', 'tutorial', x['source'])

    model = Model(os.path.join(os.getenv("HOME"), 'plom_test_model'), c, p, l)

    print model.par_fixed
    print model.par_proc

    print 'parameters: ', model.get_par_id()

    print 'order of context elements: ', model.get_ts_id(), model.get_cac_id()

    print model.ts_id
    print model.cac_id
    print model.obs_var
    print model.map_ts_obs

    model.prepare()
    model.write_settings()
    model.code()
    model.compile()

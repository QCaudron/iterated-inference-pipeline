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
import sys

class Cmodel:

    """
    parse a JSON model description:

    -remove par_fixed from par_obs and par_proc
    -replace N by sum_SV depending the context
    """

    def __init__(self, context, process, link,  **kwargs):

        self.op = set(['+', '-', '*', '/', ',', '(', ')']) ##!!!CAN'T contain square bracket '[' ']'
        self.reserved = set(['N', 'x', 'U'])
        self.special_functions = set(['terms_forcing', 'step', 'step_lin', 'sin', 'cos', 'correct_rate'])

        ###########################################################################
        self.context = copy.deepcopy(context)


        ##list [{'id':'X'}, ...] => ['X', ...]

        self.remainder = None
        self.par_sv = []
        for x in process['state']:
            if x.get('tag') == 'remainder':
                self.remainder = x['id']
            else:
                self.par_sv.append(x['id'])


        self.universes = ['U', self.remainder]

        ##par fixed !! N **HAS** to be first if it exists
        self.par_fixed = [x['id'] for x in context.get('metadata', []) if x['id'] == 'N'] #does N exist, if so it is first
        self.par_fixed += [x['id'] for x in (context.get('data', []) + context.get('metadata', [])) if x['id'] != 'data' and x['id'] != 'N'] ##all the rest

        ##remove par_fixed from par_proc and par_obs (be sure to conserve the original order so don't use set')
        self.par_proc = [x['id'] for x in process['parameter'] if x['id'] not in self.par_fixed]
        self.par_obs = [x['id'] for x in link['observation'][0]['parameter'] if x['id'] not in self.par_fixed]

        ##############################
        ## proc_model
        ##############################
        self.proc_model = copy.deepcopy(process['model'])

        ##add white_noise properties to proc_model reactions.
        self.white_noise = []
        if 'white_noise' in process:
            for j, x in enumerate(process['white_noise']):
                mynoise =  {'name': 'white_noise__' + str(j),
                            'sd': x['sd']}

                self.white_noise.append(mynoise)

                for i, r in enumerate(self.proc_model):
                    for nr in x['reaction']:
                        if r['from'] == nr['from'] and r['to'] == nr['to'] and ( ('rate' in nr and r['rate'] == nr['rate']) or ('rate' not in nr)):
                            self.proc_model[i]['white_noise'] = mynoise


        ##############################
        ## obs_model and drift
        ##############################
        self.obs_model = copy.deepcopy(link['observation'][0]['model'])

        #drift
        self.drift_par_proc = []
        self.vol_par_proc = []
        if 'diffusion' in process:
            for x in process['diffusion']:
                self.drift_par_proc.append(x['parameter'])
                self.vol_par_proc.append(x['volatility'])

        #to do support for drift in observation model
        self.drift_par_obs = []
        self.vol_par_obs = []

        #get par_fixed involved in the obs_model (par_fixed_obs):
        self.par_fixed_obs = []
        for k, v in self.obs_model.iteritems():
            if k != 'distribution' and k != 'comment':
                elements = self.change_user_input(v)
                for e in elements:
                    if e in self.par_fixed:
                        self.par_fixed_obs.append(e)

        self.par_fixed_obs = set(self.par_fixed_obs)

        self.drift_var = ['drift__par_proc__' + x for x in self.drift_par_proc]
        self.drift_var += ['drift__par_obs__' + x for x in self.drift_par_obs]


        ##############################
        ## observed variables
        ##############################

        ####IMPORTANT: we sort obs_var so that incidences are first
        sorted_obs_var= copy.deepcopy(link['observed'])
        sorted_obs_var.sort(key=lambda x: 0 if isinstance(x['definition'][0], dict) else 1)

        self.obs_var = []
        self.obs_var_def = []

        ##do a hash to quickly add rate for obs_var_def
        pm = {}
        for v in self.proc_model:
            pm[(v['from'], v['to'])] = v['rate']

        for x in sorted_obs_var:
            self.obs_var.append(x['id'])
            ##add rates if not present to incidences of obs_var_def
            mydef = copy.deepcopy(x['definition'])
            for d in mydef:
                if isinstance(d, dict) and 'rate' not in d:
                    d['rate'] = pm[(d['from'], d['to'])]

            self.obs_var_def.append(mydef)


        ##add white_noise properties to incidences reactions of obs_var_def.
        if 'white_noise' in process:
            for i, inc_prev_list in enumerate(self.obs_var_def):
                for k, r in enumerate(inc_prev_list):
                    if isinstance(r, dict): #only for incidence
                        for j, x in enumerate(process['white_noise']):
                            for nr in x['reaction']:
                                if r['from'] == nr['from'] and r['to'] == nr['to'] and (('rate' in nr and r['rate'] == nr['rate']) or ('rate' not in nr)):
                                    self.obs_var_def[i][k]['white_noise'] = {'name': 'white_noise__' + str(j),
                                                                             'sd': x['sd']}

                                    
        ##We treat reaction starting from remainder as reaction starting from U that is rate -> rate * from size. It results in simpler code in Ccoder.py
        ##if no remainder, we replace 'N' by sum of par_sv
        ##else, we replace remainder by N - sum(par_sv) in the rates (and in the rates ONLY)

        resolve_N = lambda x: '({0})'.format('+'.join(self.par_sv)) if x == 'N' else x
        resolve_remainder = lambda x: '(N-{0})'.format('-'.join(self.par_sv)) if x == self.remainder else x

        if not self.remainder:
            for k, v in self.obs_model.iteritems():
                if k != 'distribution':
                    self.obs_model[k] = ''.join(map(resolve_N, self.change_user_input(v)))

        for i, m in enumerate(self.proc_model):
            if self.remainder:
                if self.proc_model[i]['from'] == self.remainder:
                    self.proc_model[i]['rate'] = '({0})*{1}'.format(self.proc_model[i]['rate'], self.remainder)

                self.proc_model[i]['rate'] = ''.join(map(resolve_remainder, self.change_user_input(m['rate'])))

            else:
                self.proc_model[i]['rate'] = ''.join(map(resolve_N, self.change_user_input(m['rate'])))

        for x in self.obs_var_def:
            for d in x:
                if isinstance(d, dict): #incidence
                    if self.remainder:
                        if d['from'] == self.remainder:
                            d['rate'] = '({0})*{1}'.format(d['rate'], self.remainder)

                        d['rate'] = ''.join(map(resolve_remainder, self.change_user_input(d['rate'])))

                    else:
                        d['rate'] = ''.join(map(resolve_N, self.change_user_input(d['rate'])))                        

        ##TODO other models (spate or age structure)



    def change_user_input(self, reaction):
        """transform the reaction in smtg that we can parse in a programming language:
        example: change_user_input('r0*2*correct_rate(v)') -> ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')']"""

        myreaction=reaction.replace(' ','') ##get rid of whitespaces
        mylist=[]
        mystring=''

        for i in range(len(myreaction)):

            if myreaction[i] in self.op :
                if len(mystring)>0:
                    mylist.append(mystring)
                    mystring=''
                mylist.append(myreaction[i])
            else:
                mystring += myreaction[i]

        if len(mystring)>0: ##the string doesn't end with an operator
            mylist.append(mystring)

        return mylist






if __name__=="__main__":


    print 'Cf test_Cmodel.py for tests'

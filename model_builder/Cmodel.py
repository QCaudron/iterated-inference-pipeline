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
        self.reserved = set(['sum_SV', 'N', 'prop', 'x', 'U'])
        self.special_functions = set(['terms_forcing', 'step', 'step_lin', 'sin', 'cos', 'correct_rate'])

        ###########################################################################
        self.context = copy.deepcopy(context)


        ##list [{'id':'X'}, ...] => ['X', ...]

        self.remainder = {}
        remainder_key = None
        self.par_sv = []
        for x in process['state']:
            if x.get('tag') == 'remainder':
                remainder_key = x['id']
            else:
                self.par_sv.append(x['id'])

        if remainder_key:
            self.remainder[remainder_key] = '(N -' + '-'.join(self.par_sv) + ')'

        self.universes = ['U'] + self.remainder.keys()

        ##par fixed !! N **HAS** to be first if it exists
        self.par_fixed = [x['id'] for x in context.get('metadata', []) if x['id'] == 'N'] #does N exist, if so it is first
        self.par_fixed += [x['id'] for x in (context.get('data', []) + context.get('metadata', [])) if x['id'] != 'data' and x['id'] != 'prop' and x['id'] != 'N'] ##all the rest

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


                                    
        ##resolve remainder (e.g R -> (N-S-I))        
        ## treat reaction starting from remainder as reaction starting from U that is rate -? rate * from size (simpler code in Ccoder)
        ## e.g, R is remainder and  R -> S, g | R - > s, g*(N-S-I)

        ##OR resolve the population size: (replace 'N' by 'sum_SV' if no remainder)
        if self.remainder:

            resolve_remainder = lambda x: self.remainder[x] if x in self.remainder else x
            resolve_N = lambda x: 'sum_SV' if x == 'N' else x

            for k, v in self.obs_model.iteritems():
                if k != 'distribution':
                    self.obs_model[k] = ''.join(map(resolve_N, self.change_user_input(v)))

            for i, m in enumerate(self.proc_model):
                if self.remainder:
                    self.proc_model[i]['rate'] = ''.join(map(resolve_remainder, self.change_user_input(m['rate'])))
                    if self.proc_model[i]['from'] in self.remainder:
                        self.proc_model[i]['rate'] = '({0})*{1}'.format(self.proc_model[i]['rate'], self.remainder[ self.proc_model[i]['from'] ])

                else:
                    self.proc_model[i]['rate'] = ''.join(map(resolve_N, self.change_user_input(m['rate'])))

            for x in self.obs_var_def:
                for d in x:
                    if isinstance(d, dict):
                        if self.remainder:
                            d['rate'] = ''.join(map(resolve_remainder, self.change_user_input(d['rate'])))
                            if d['from'] in self.remainder:
                                d['rate'] = '({0})*{1}'.format(d['rate'], self.remainder[ d['from'] ])

                        else:
                            d['rate'] = ''.join(map(resolve_N, self.change_user_input(d['rate'])))



        ##TODO other models (spate or age structure)



    def change_user_input(self, reaction):
        """transform the reaction in smtg that we can parse in a programming language:
        example: change_user_input('r0*2*correct_rate(v)') -> ['r0', '*', '2', 'correct_rate', '(', 'v', ')']"""

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

    """
    test substitution of N in either sum_SV or N
    """

    m = {}
    m['state'] = [{'id':'S'}, {'id':'I'}, {'id': 'R', 'tag': 'remainder'}]
    m['parameter'] = [{'id':'r0'}, {'id':'v'}, {'id':'l'}, {'id':'e'}, {'id':'d'}, {'id':'sto'}, {'id':'alpha'}, {'id':'mu_b'}, {'id':'mu_d'}, {'id':'vol'}, {'id':'g'}]


    m['model'] = [ {'from': 'U', 'to': 'S',  'rate': 'mu_b*N'},
                   {'from': 'S', 'to': 'E',  'rate': 'r0/N*v*(1.0+e*sin_t(d))*I', "tag": 'transmission'},

                   {'from': 'E', 'to': 'I', 'rate': '(1-alpha)*correct_rate(l)'},
                    ##Here we split the reaction from E->U as we only observe a subpart
                   {'from': 'E', 'to': 'U',  'rate': 'alpha*correct_rate(l)'},
                   {'from': 'E', 'to': 'U',  'rate': 'mu_d'},

                   {'from': 'S', 'to': 'U',  'rate': 'mu_d'},
                   {'from': 'R', 'to': 'S',  'rate': 'g'},
                   {'from': 'I', 'to': 'R', 'rate': '(1-alpha)*correct_rate(v)'},
                   {'from': 'I', 'to': 'U',  'rate': 'alpha*correct_rate(v) + mu_d'} ]

    m['diffusion'] = [{'parameter':'r0',
                       'volatility': 'vol',
                       'drift': 0.0}]

    m['white_noise'] = [{'reaction': [{'from':'S', 'to': 'E'}],
                         'sd': 'sto'}]

    ##context elements needed for Cmodel
    c = {}

    c['data'] = [{'id': 'data'}, {'id': 'prop'}]

    c['metadata'] = [{'id': 'mu_b'}, {'id': 'mu_d'}, {'id': 'N'}]

    ##link elements needed for Cmodel
    l = {}
    l['observed'] =  [{"id": "Prev",     "definition": ["I"], "model_id": "common"},
                      {"id": "SI",       "definition": ["S", "I"], "model_id": "common"},
                      ##we have to specify a rate to the incidence E->U as we only observed a subpart of this reaction
                      {"id": "Inc_out",  "definition": [{"from":"I", "to":"R"}, {"from":"E", "to":"U", 'rate': "mu_d"}], "model_id": "common"},
                      {"id": "Inc_in",   "definition": [{"from":"S", "to":"E"}], "model_id": "common"},
                      {"id": "Inc_x",   "definition": [{"from":"R", "to":"S"}], "model_id": "common"}]


    l["observation"] = [{"id": "common", 
                         "parameter": [{"id": "rep","comment": "reporting rate"}, 
                                       {"id": "phi",  "comment": "over-dispertion"}],
                         "model": {"distribution": "discretized_normal",
                                   "mean": "rep*prop*x",
                                   "var": "rep*(1.0-rep)*prop*x + (rep*phi*prop*x)**2"}}]

    test_model = Cmodel(c, m, l)

    print "par_fixed"
    print test_model.par_fixed

    print "par_proc"
    print test_model.par_proc

    print "state variables"
    print test_model.par_sv

    print "drift var"
    print test_model.drift_par_proc
    print test_model.vol_par_proc
    print test_model.drift_par_obs
    print test_model.vol_par_obs
    print test_model.drift_var

    print "\nprocess model"
    tmp_print = test_model.proc_model
    for line in tmp_print:
        wn = line['white_noise'] if 'white_noise' in line else ''
        print line['from'], line['to'], line['rate'], wn

    print "\nobserved variable definition"
    tmp_print = test_model.obs_var_def
    for line in tmp_print:
        print line

    print "\ntest_model.obs_model"
    print test_model.obs_model

    print "\ntest_remainder"
    print test_model.remainder

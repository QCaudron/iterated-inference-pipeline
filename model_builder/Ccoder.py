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
from Cmodel import Cmodel
from sympy import diff, Symbol, sympify, simplify
from sympy.printing import ccode
import copy

class Ccoder(Cmodel):
    """write the C code from the user input coming from the web interface..."""

    def __init__(self, context, process, link,  **kwargs):
        Cmodel.__init__(self, context, process, link,  **kwargs)

        self.set_par_fixed = set(self.par_fixed)
        self.all_par = self.par_sv + [self.remainder] + self.drift_par_proc + self.par_proc +  self.par_obs + ['x'] + ['N'] + ['t']

    def toC(self, term, no_correct_rate):

        if term in self.par_sv:
            return 'X[ORDER_{0}*N_CAC+cac]'.format(term)

        elif term in self.set_par_fixed:
            if term in self.par_fixed_obs:
                return 'covar[ORDER_{0}][nn][ts]'.format(term)
            else:
                return 'covar[ORDER_{0}][nn][cac]'.format(term)

        elif term in self.par_proc :
            if term in self.drift_par_proc:
                return 'drifted[ORDER_drift__par_proc__{0}][cac]'.format(term)
            else:
                return 'par[ORDER_{0}][routers[ORDER_{0}]->map[cac]]'.format(term)

        elif term in self.par_obs:
            ##TODO support for drift for parameters of the observation process
            return 'par[ORDER_{0}][routers[ORDER_{0}]->map[ts]]'.format(term)

        elif term == 'correct_rate' and no_correct_rate:
            return ''

        else: ##r is an operator or x
            return term


    def generator_C(self, term, no_correct_rate):

        terms = self.change_user_input(term)

        ind = 0
        Cterm = ''
        while (ind < len(terms)):

            if terms[ind] in self.special_functions:
                myf = terms[ind]

                Cterm += self.toC(myf, no_correct_rate) + '('
                ind += 2 #skip first parenthesis

                pos = 1 #counter for open parenthesis
                while pos > 0:
                    if terms[ind] == '(':
                        pos += 1
                    if terms[ind] == ')':
                        pos -= 1

                    if pos >0:
                        Cterm += self.toC(terms[ind], no_correct_rate)
                        ind += 1

                ##add extra terms (no whitespace)
                if myf == 'terms_forcing':
                    Cterm += ',t,p_data,cac'
                elif myf == 'step':
                    Cterm += ',t'
                elif myf == 'step_lin':
                    Cterm += ',t'
                elif myf == 'correct_rate' and not no_correct_rate:
                    Cterm += ',dt'

                ##close bracket
                Cterm += terms[ind]
                ind += 1

            else:
                Cterm += self.toC(terms[ind], no_correct_rate)
                ind += 1

        return Cterm


    def make_C_term(self, term, no_correct_rate, derivate=None, human=False):

        """transform a term into its plom C expression OR the
        plom C expression of its derivate, differentiating
        against the derivate (if derivate not None)

        """

        #prefix all the state variable and parameters by plom___ to
        #avoid namespace collision with Sympy as QCOSINE letters are
        #used by SymPy

        myterm = self.change_user_input(term)
        safe = ''

        for r in myterm:
            if r in self.all_par:
                safe += 'plom___' + r
            else:
                safe += r

        if derivate:
            sy = Symbol(str('plom___'+derivate))
            pterm = diff(sympify(safe), sy)
        else:
            pterm = sympify(safe)

        #remove the plom___ prefix                
        #term = ccode(simplify(pterm)).replace('plom___', '') ##NOTE simplify is just too slow to be used...
        term = ccode(pterm).replace('plom___', '')

        #make the plom C expression
        if human:
            return term
        else:            
            return self.generator_C(term, no_correct_rate)


    def print_obs_prev(self):
        """generate C code to track the observed **prevalence** """

        Clist = []

        for i in range(len(self.obs_var_def)):
            if not isinstance(self.obs_var_def[i][0], dict): ##prevalence                
                prev = '+'.join(map(lambda x: '(N-{0})'.format('-'.join(self.par_sv)) if x == self.remainder else x, self.obs_var_def[i]))
                Clist.append(self.make_C_term(prev, False))

        return Clist


    def obs_inc_step_psr(self):
        """generate C code to compute the dynamic of the observed
        **incidence** in case of stochastic models (euler multinomial)
        and put in into

        Clist = [{'true_ind_obs':'i (i count both incidence and prevalence (N_OBS_ALL))',
                  'right_hand_side':'f(X[i,c,ac]'}]

        """
        Clist = []

        for i in range(len(self.obs_var_def)):
            true_ind_obs = str(i)
            right_hand_side=''

            if isinstance(self.obs_var_def[i][0], dict): ##incidence

                for j in range(len(self.obs_var_def[i])):
                    id_out = [self.proc_model.index(r) for r in self.proc_model if ((r['from'] == self.obs_var_def[i][j]['from']) and (r['to'] == self.obs_var_def[i][j]['to']) and (r['rate'] == self.obs_var_def[i][j]['rate']))]
                    for o in id_out:
                        myexit = [r for r in self.proc_model if r['from']==self.proc_model[o]['from']]
                        right_hand_side += ' + p_calc->inc[ORDER_{0}][cac][{1}]'.format(self.obs_var_def[i][j]['from'], myexit.index(self.proc_model[o]))

                Clist.append({'true_ind_obs':true_ind_obs, 'right_hand_side':right_hand_side})

        return Clist


    def print_build_psr(self):
        Clist = []
        for s in self.par_sv + ['U', self.remainder]:
            nbreac = len([r for r in self.proc_model if r['from']==s]) +1 ##+1 to stay in the same compartment or to declare smtg in case of no reaction (not super clean but makes C code easier...)
            Clist.append({'state':s, 'nb_reaction': nbreac})

        return Clist


    def step_psr(self):

        """
        prob and update for Poisson with stochastic rate step function

        prob general structure:

        sum=...;
        if(sum>0.0){
            prob[0]=(1-exp(-sum))*(rate/sum);
            ...
            prob[last]=1-sum(prob);
        }
        else{
            prob[0]=0;
            ...
            prob[last]=1;
        }
        we need the if statement to avoid division by 0
        """

        ###########
        ## prob  ##
        ###########
        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        ##make the rates noisy (if needed e.g r
        rates = set()
        for r in proc_model:
            if r['from'] not in ['U', self.remainder]:
                myrate = r['rate']
                if 'white_noise' in r:
                    myrate = '({0})*{1}'.format(myrate, r['white_noise']['name'])

                rates.add(myrate)

        rates = list(rates)
        caches = map(lambda x: self.make_C_term(x, False), rates)
        sf = self.cache_special_function_C(caches)

        for r in proc_model:
            if r['from'] not in ['U', self.remainder]:
                myrate = r['rate']
                if 'white_noise' in r:
                    myrate = '({0})*{1}'.format(myrate, r['white_noise']['name'])

                r['ind_cache'] = rates.index(myrate)

        Ccode=''

        for s in self.par_sv:
            myexit = [r for r in proc_model if r['from'] == s]
            exitlist=[]

            if len(myexit)>0:

                for e in myexit:
                    exitlist.append('_r[{0}]*dt'.format(e['ind_cache']))

                Csum= 'sum = ' + '+'.join(exitlist) + ';\n' #.join() doesn't add '+' if len(list)==1
                Ccode += Csum+ 'if(sum>0.0){\none_minus_exp_sum = (1.0-exp(-sum));\n'
                Cprob=''
                sumprob='1.0'
                for reacnb in range(len(exitlist)):
                    Cprob += 'p_calc->prob[ORDER_{0}][{1}] = one_minus_exp_sum*(({2})/sum);\n'.format(s, reacnb, exitlist[reacnb])
                    sumprob += ' - p_calc->prob[ORDER_{0}][{1}]'.format(s, reacnb)

                Cprob += 'p_calc->prob[ORDER_{0}][{1}] = '.format(s,len(exitlist)) + sumprob + ';\n'
                Ccode += Cprob+ '}\n'
                Ccode +='else{\n'

                Celse=''
                for reacnb in range(len(exitlist)):
                    Celse += 'p_calc->prob[ORDER_{0}][{1}] = 0.0;\n'.format(s, reacnb)

                Celse += 'p_calc->prob[ORDER_{0}][{1}] = 1.0;\n'.format(s,len(exitlist))+'}\n\n'

                Ccode += Celse

        ############
        ## update ##
        ############

        incDict = dict([(x,'') for x in self.par_sv])

        for s in self.par_sv: ##stay in the same compartment
            myexit = [r for r in self.proc_model if r['from'] == s]
            if len(myexit)>0: ##only if you can exit from this compartment in this case the remaining has a sense
                incDict[s] += 'p_calc->inc[ORDER_{0}][cac][{1}]'.format(s, len(myexit))
            else:
                incDict[s] += 'X[ORDER_{0}*N_CAC+cac]'.format(s)

        for s in self.par_sv: #come in from other compartments
            myinput = [r for r in self.proc_model if r['from'] == s]
            for nbreac in range(len(myinput)):
                if myinput[nbreac]['to'] not in ['U', self.remainder]: ##we exclude deaths or transitions to remainder in the update
                    incDict[myinput[nbreac]['to']] += ' + p_calc->inc[ORDER_{0}][cac][{1}]'.format(myinput[nbreac]['from'], nbreac)


        ##we add flow from ['U', self.remainder] (Poisson term). We want to cache those flow so that the incidences can be computed
        poisson = []
        for s in ['U', self.remainder]:
            reac_from_univ = [r for r in self.proc_model if r['from'] == s]
            for nbreac in range(len(reac_from_univ)):
                myrate = self.make_C_term(reac_from_univ[nbreac]['rate'], False)
                if 'white_noise' in reac_from_univ[nbreac]:
                    myrate = '({0})*{1}'.format(myrate, reac_from_univ[nbreac]['white_noise']['name'])

                poisson.append('p_calc->inc[ORDER_{0}][cac][{1}] = gsl_ran_poisson(p_calc->randgsl, ({2})*dt)'.format(s, nbreac, myrate))
                incDict[reac_from_univ[nbreac]['to']] += ' + p_calc->inc[ORDER_{0}][cac][{1}]'.format(s, nbreac)

        Cstring=''
        for s in self.par_sv:
            Cstring += 'X[ORDER_{0}*N_CAC+cac] = {1};\n'.format(s, incDict[s])


        return {'code': Ccode, 'caches': caches, 'sf': sf, 'poisson': poisson, 'update': Cstring}


    def psr_multinomial(self):
        draw = []
        for s in self.par_sv:
            nbexit = len([r for r in self.proc_model if r['from']==s])
            if nbexit>0:
                draw.append({'state': s, 'nb_exit': nbexit+1}) ##+1 to stay in the compartment

        return draw


    def cache_special_function_C(self, caches_C, sf=None, prefix='_sf'):
        """caches_C: List of cached expression in C
        caches_C is modified in place
        sf: an optional list of unique special function to be cached
        returns sf (created if sf input is None)
        """

        if not sf:
            sf = []
            for term in caches_C:
                if any([x in term for x in self.special_functions]):
                    terms = self.change_user_input(term)
                    ind = 0
                    while (ind < len(terms)):
                        if terms[ind] in self.special_functions:
                            f = terms[ind] + '('
                            ind += 2 #skip first parenthesis
                            pos = 1 #counter for open parenthesis
                            while pos > 0:
                                if terms[ind] == '(':
                                    pos += 1
                                if terms[ind] == ')':
                                    pos -= 1

                                f += terms[ind]
                                ind +=1
                        
                            sf.append(f)
                        else:
                            ind += 1

            sf = list(set(sf))

        for i, term in enumerate(caches_C):
            if any([x in term for x in self.special_functions]):
                for s in sf:
                    caches_C[i] = caches_C[i].replace(s, prefix + '[{0}]'.format(sf.index(s)))

        return sf


    def step_ode_sde(self):
        """
        Generates ODE and SDEs
        note: sf are used in self.jac() for Lyapunov exp computations
        """

        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        odeDict = dict([(x, []) for x in self.par_sv])

        rates = list(set(r['rate'] for r in proc_model))

        caches = map(lambda x: self.make_C_term(x, True), rates)
        sf = self.cache_special_function_C(caches, prefix='_sf[cac]')

        for i, r in enumerate(proc_model):
            r['ind_cache'] = rates.index(r['rate'])
            r['ind_dem_sto'] = i


        def get_rhs_term(sign, cached, reaction):
            if 'white_noise' in reaction:
                noise_name = reaction['white_noise']['name']
                noise_sd = self.toC(reaction['white_noise']['sd'], False)
            else:
                noise_name = None
                noise_sd= None

            return {'sign': sign, 'term': cached, 'noise_name': noise_name, 'noise_sd': noise_sd, 'ind_dem_sto': reaction['ind_dem_sto']}


        ################################
        ##Dynamic of the state variables
        ################################

        ##outputs
        for r in proc_model:
            if r['from'] not in ['U', self.remainder]:
                cached = '_r[cac][{0}]*X[ORDER_{1}{2}]'.format(r['ind_cache'], r['from'], '*N_CAC+cac')
                odeDict[r['from']].append(get_rhs_term('-', cached, r))

        ##inputs
        for r in proc_model:
            if r['to'] not in ['U', self.remainder]:
                if r['from'] not in ['U', self.remainder]:
                    cached = '_r[cac][{0}]*X[ORDER_{1}{2}]'.format(r['ind_cache'], r['from'], '*N_CAC+cac')
                else:
                    cached= '_r[cac][{0}]'.format(r['ind_cache'])

                odeDict[r['to']].append(get_rhs_term('+', cached, r))


        #######################################
        ##Dynamic of the observed **incidence**
        #######################################

        obs_list = []

        for i in range(len(self.obs_var_def)):
            eq = []

            if isinstance(self.obs_var_def[i][0], dict): ##incidence
                for j in range(len(self.obs_var_def[i])):
                    id_out = [proc_model.index(r) for r in proc_model if ((r['from'] == self.obs_var_def[i][j]['from']) and (r['to'] == self.obs_var_def[i][j]['to']) and (r['rate'] == self.obs_var_def[i][j]['rate'])) ]
                    for o in id_out:
                        reaction = proc_model[o]
                        if self.obs_var_def[i][j]['from'] in ['U', self.remainder]:
                            cached = '_r[cac][{0}]'.format(reaction['ind_cache'])
                        else:
                            cached = '_r[cac][{0}]*X[ORDER_{1}{2}]'.format(reaction['ind_cache'], self.obs_var_def[i][j]['from'], '*N_CAC+cac')

                        eq.append(get_rhs_term('+', cached, reaction))

                obs_list.append({'index':i, 'eq': eq})



        ##############################################################################################################
        ##we create the ODE and  4 versions of the SDE system (no_dem_sto, no_env_sto, no_dem_sto_no_env_sto and full)
        ##############################################################################################################
        unique_noises_names = [x['name'] for x in self.white_noise]
        dem_sto_names = ['dem_sto__' +str(i) for i, x in enumerate(self.proc_model)]

        def eq_dem_env(eq_list):
            eq = ''  #deter skeleton
            dem = '' #demographic stochasticity
            env = '' #env stochasticity

            for x in eq_list:
                eq += ' {0} ({1})'.format(x['sign'], x['term'])

                #dem sto
                dem += '{0} sqrt(({1}))*dem_sto__{2}[cac]'.format(x['sign'], x['term'], x['ind_dem_sto'])

                #env sto
                if x['noise_name']:
                    env += '{0} ({1})*{2}*{3}[cac]'.format(x['sign'], x['term'], x['noise_sd'], x['noise_name'])

            return (eq, dem, env)


        func = {'no_dem_sto': {'proc': {'system':[], 'noises': unique_noises_names},
                               'obs': []},
                'no_env_sto': {'proc': {'system':[], 'noises': dem_sto_names},
                               'obs': []},
                'full': {'proc': {'system':[], 'noises': dem_sto_names + unique_noises_names},
                         'obs': []},
                'no_dem_sto_no_env_sto': {'proc':{'system':[], 'noises':[]},
                                          'obs':[]},
                'ode': {'proc':{'system':[], 'noises':[]},
                        'obs':[]}}


        #state variables
        for i, s in enumerate(self.par_sv):

            eq, dem, env = eq_dem_env(odeDict[s])
            if env:
                env = '+ ' + env

            #TODO get rid of the 'dt' for Euler Maruyama (should be handled on the C side as it is the case for sqrt(dt) for the stochastic part)'
            func['ode']['proc']['system'].append({'index': i, 'eq': eq})
            func['no_dem_sto_no_env_sto']['proc']['system'].append({'index': i, 'eq': '({0})*dt'.format(eq)})
            func['no_dem_sto']['proc']['system'].append({'index': i, 'eq': '({0})*dt {1}'.format(eq, env)})
            func['no_env_sto']['proc']['system'].append({'index': i, 'eq': '({0})*dt + {1}'.format(eq, dem)})
            func['full']['proc']['system'].append({'index': i, 'eq': '({0})*dt + {1} {2}'.format(eq, dem, env)})

        #observed incidence
        for myobs in obs_list:

            eq, dem, env = eq_dem_env(myobs['eq'])
            if env:
                env = ' + ' + env

            #TODO get rid of the 'dt' for Euler Maruyama (should be handled on the C side as it is the case for sqrt(dt) for the stochastic part)'
            func['ode']['obs'].append({'index': myobs['index'], 'eq': eq})
            func['no_dem_sto_no_env_sto']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt'.format(eq)})
            func['no_dem_sto']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt {1}'.format(eq, env)})
            func['no_env_sto']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt {1}'.format(eq, dem)})
            func['full']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt + {1} {2}'.format(eq, dem, env)})


        return {'func': func, 'caches': caches, 'sf': sf}




    def print_order(self):
        order_univ = []
        N_PAR_SV = len(self.par_sv)
        for i, X in enumerate(['U', self.remainder]):
            order_univ.append({'name': X, 'order': N_PAR_SV+i})

        return {'var': self.par_sv + self.par_proc + self.par_obs, 'drift': self.drift_var, 'data': self.par_fixed, 'universe': order_univ}


    def print_like(self):

        ##WARNING right now only the discretized normal is supported.
        ##TODO: generalization for different distribution

        return {'mean':self.make_C_term(self.obs_model['mean'], False),
                'var':self.make_C_term(self.obs_model['var'], False)}



    def jac(self, sf_jac_only):
        """compute jacobian matrix of the process model (including
        observed variable) using Sympy


        sf_jac_only: list of cached special function generated by
        self.print_ode() used to get the index of caches_C for the
        jacobian matrix of simulation methods

        """

        my_model = copy.deepcopy(self.proc_model)
        odeDict = dict([(x,'') for x in self.par_sv])

        ##outputs
        for r in my_model:
            if r['from'] not in ['U', self.remainder]:
                rate= ' - (({0})*{1})'.format(r['rate'], r['from'])
                odeDict[r['from']] += rate

        ##inputs
        for r in my_model:
            if r['to'] not in ['U', self.remainder]:
                if r['from'] not in ['U', self.remainder]:
                    rate= ' + (({0})*{1})'.format(r['rate'], r['from'])
                    odeDict[r['to']] += rate
                else:
                    rate= ' + ({0})'.format(r['rate'])
                    odeDict[r['to']] += rate


        ##observed equations
        obsList = []

        for i in range(len(self.obs_var_def)):
            eq = ''

            if not isinstance(self.obs_var_def[i][0], dict): ##prevalence: we potentialy sum the prevalence eq. stored in odeDict
                for j in range(len(self.obs_var_def[i])):
                    if self.obs_var_def[i][j] == self.remainder:
                        for s in self.par_sv:
                            eq += '- ( ' + odeDict[s] + ' )'
                    else:
                        eq += odeDict[self.obs_var_def[i][j]]

            else: ##incidence
                for j in range(len(self.obs_var_def[i])):
                    id_out = [self.proc_model.index(r) for r in self.proc_model if ((r['from'] == self.obs_var_def[i][j]['from']) and (r['to'] == self.obs_var_def[i][j]['to']) and (r['rate'] == self.obs_var_def[i][j]['rate'])) ]
                    for o in id_out:
                        reaction = my_model[o]
                        if self.obs_var_def[i][j]['from'] in ['U', self.remainder]:
                            eq += ' + ({0})'.format(reaction['rate'])
                        else:
                            eq += ' + (({0})*{1})'.format(reaction['rate'], self.obs_var_def[i][j]['from'])

            obsList.append(eq)


        ####################
        ### Jacobian
        ####################

        ##derive process model equations (odeDict) per par_sv
        caches = []
        caches_jac_only = []

        jac = []
        jac_only = []
        jac_drift = []

        for s in range(len(self.par_sv)):
            jac.append([])
            jac_only.append([])

            if self.drift_var:
                jac_drift.append([])

            for sy in self.par_sv:
                Cterm = self.make_C_term(odeDict[self.par_sv[s]], True, derivate=sy)
                jac[s].append(Cterm)
                jac_only[s].append(Cterm)
                caches.append(Cterm)
                caches_jac_only.append(Cterm)

            #see doc of kalman.c drift_derivative()
            for sy in self.drift_par_proc:
                Cterm = self.make_C_term(odeDict[self.par_sv[s]], True, derivate=sy)
                jac_drift[s].append({'value': Cterm,
                                     'der': self.make_C_term(sy, True),
                                     'name': sy})
                caches.append(Cterm)

        ##derive observation equations (obsList) per par_sv
        jac_obs = []
        jac_obs_drift = []

        for o in range(len(obsList)):
            jac_obs.append([])
            if self.drift_var:
                jac_obs_drift.append([])

            for sy in self.par_sv:
                Cterm = self.make_C_term(obsList[o], True, derivate=sy)

                jac_obs[o].append(Cterm)
                caches.append(Cterm)

            #see doc of kalman.c drift_derivative()
            for sy in self.drift_par_proc:
                Cterm = self.make_C_term(obsList[o], True, derivate=sy)
                jac_obs_drift[o].append({'value': Cterm,
                                         'der': self.make_C_term(sy, True),
                                         'name': sy})
                caches.append(Cterm)


        ##cache rates and remove duplicates
        caches = list(set(caches))
        caches_jac_only = list(set(caches_jac_only))

        ##replace with index of caches (will be _r[index] in C)
        for s in range(len(self.par_sv)):
            for i in range(len(self.par_sv)):
                Cterm = jac[s][i]
                jac[s][i] = caches.index(Cterm)
                jac_only[s][i] = caches_jac_only.index(Cterm)

            for i in range(len(self.drift_par_proc)):
                jac_drift[s][i]['value'] = caches.index(jac_drift[s][i]['value'])


        for o in range(len(obsList)):
            for i in range(len(self.par_sv)):
                jac_obs[o][i] = caches.index(jac_obs[o][i])

            for i in range(len(self.drift_par_proc)):
                jac_obs_drift[o][i]['value'] = caches.index(jac_obs_drift[o][i]['value'])


        ##special function that have to be cached (caches is transformed by self.cache_special_function_)
        sf = self.cache_special_function_C(caches, prefix='_sf[cac]')
        ##for jac_only (used for Lyapunov exp computations only, sf is shared with the one of print_ode. We just update caches_jac_only)
        self.cache_special_function_C(caches_jac_only, sf=sf_jac_only, prefix='_sf[cac]')

        return {'jac_only':jac_only,
                'jac':jac,
                'jac_obs':jac_obs,
                'jac_drift':jac_drift,
                'jac_obs_drift':jac_obs_drift,
                'caches': caches,
                'sf': sf,
                'caches_jac_only': caches_jac_only}


    def der_mean_proc_obs(self):
        """compute jacobian matrix of the mean of the obs process (assumed to be Gaussian) using Sympy"""

        return self.make_C_term(self.obs_model['mean'], True, derivate='x')

    def der2_mean_proc_obs(self):
        """compute the second derivative of the mean of the obs process (assumed to be Gaussian) using Sympy"""

        first = self.make_C_term(self.obs_model['mean'], True, derivate='x', human=True)

        return self.make_C_term(first, True, derivate='x')


    def eval_Q(self, debug = False):
        """

        The construction of Qsv is based on three levels:
         - states: state variables and observations (s)
         - reactions (r)
         - noise terms (n)

        At the reaction level, Qr is a two-blocks diagonal matrix: Qr_dem and Qr_env.
        Qr_dem corresponds to demographic noise and has reaction rates on the diagonal.
        Qr_env corresponds to white noises. It is built from Qn through Lr.
        Qn is a diagonal matrix which diagonal terms correspond to squarred amplitude of white noises.
        The stoechiometric matrices L are used to switch from one level to another:
              Qr_env = Lr Qn Lr'  and Qs = Ls Qr Ls'
        
        In particular, Lr has reaction rates in term (i,j) if reaction i is concerned by white noise j.
        Ls has +1 or -1 in term (i,j) if reaction j goes to or leaves from state i, and O's everywhere else.

        Note: we assume only one environmental noise term per reaction


        """
        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        N_REAC = len(proc_model)
        N_PAR_SV = len(self.par_sv)
        N_OBS = len(self.obs_var_def)

        unique_noises_names = [x['name'] for x in self.white_noise]
        N_ENV_STO_UNIQUE = len(unique_noises_names)
        
        ##add sd and order properties to noisy reactions
        N_ENV_STO = 0
        for r in proc_model:
            if 'white_noise' in r:
                r['order_env_sto_unique'] = unique_noises_names.index(r['white_noise']['name'])
                r['order_env_sto'] = N_ENV_STO
                N_ENV_STO += 1


        s = N_REAC + N_ENV_STO ##for demographic stochasticity, one independent noise term per reaction

        Ls = [[0]*s for x in range(N_PAR_SV + N_OBS)]
        Qs = [[0]*(N_PAR_SV + N_OBS) for x in range(N_PAR_SV + N_OBS)]
        Qr = [[0]*s for x in range(s)]
        Qr_dem = [[0]*s for x in range(N_REAC)]
        Qr_sto = [[0]*s for x in range(N_ENV_STO)]
        Lr = [[0]*N_ENV_STO_UNIQUE for x in range(N_ENV_STO)]
        Qn = [[0]*N_ENV_STO_UNIQUE for x in range(N_ENV_STO_UNIQUE)]


        ###########################################
        #    Create Ls and Qr_dem                 #
        ###########################################

        #state variables
        for B_dem_ind, r in enumerate(proc_model):
            is_noise = 'white_noise' in r
            if is_noise:
                B_sto_ind = N_REAC + r['order_env_sto']

            if r['from'] not in ['U', self.remainder]:
                i = self.par_sv.index(r['from'])
                Ls[i][B_dem_ind] -= 1 ##demographic stochasticity
                if is_noise:
                    Ls[i][B_sto_ind] -= 1 ##env stochasticity

                Qc_term = '({0})*{1}'.format(r['rate'], r['from'])
            else:
                Qc_term = r['rate']

            if r['to'] not in ['U', self.remainder]:
                i = self.par_sv.index(r['to'])
                Ls[i][B_dem_ind] += 1
                if is_noise:
                    Ls[i][B_sto_ind] += 1

            Qr_dem[B_dem_ind][B_dem_ind] =  Qc_term

        # observed variables
        for i in range(len(self.obs_var_def)): #(for every obs variable)

            for B_dem_ind, r in enumerate(proc_model):
                is_noise = 'white_noise' in r
                if is_noise:
                    B_sto_ind = N_REAC + r['order_env_sto']

                ##prevalence: (TODO: get rid of prevalence)
                if not isinstance(self.obs_var_def[i][0], dict):

                    # if it involves prevalence as input
                    if r['from'] in self.obs_var_def[i]:
                        Ls[N_PAR_SV + i][B_dem_ind] -= 1
                        if is_noise:
                            Ls[N_PAR_SV + i][B_sto_ind] -= 1

                    # if it involves prevalence as output
                    if r['to'] in self.obs_var_def[i]:
                        Ls[N_PAR_SV + i][B_dem_ind] += 1
                        if is_noise:
                            Ls[N_PAR_SV + i][B_sto_ind] += 1

                #incidence:
                else:
                    # for every incidence
                    for inc in self.obs_var_def[i]:
                        # if it involves incidence
                        if (r['from'] == inc['from']) and (r['to'] == inc['to']) and (r['rate'] == inc['rate']):
                            Ls[N_PAR_SV + i][B_dem_ind] += 1
                            if is_noise:
                                Ls[N_PAR_SV + i][B_sto_ind] += 1




        ############################
        ## Create Qr_env = Lr Qn Lr'
        ############################
        for r in proc_model:
            if 'white_noise' in r:
                if r['from'] not in ['U', self.remainder]:
                    Qn_term = '({0})*{1}'.format(r['rate'], r['from'])
                else:
                    Qn_term = r['rate']

                Lr[r['order_env_sto']][r['order_env_sto_unique']] = Qn_term
                Qn[r['order_env_sto_unique']][r['order_env_sto_unique']] = '({0})**2'.format(r['white_noise']['sd'])



        def matrix_product(A, B):
            if not A or not B:
                return []

            res = [[0]*len(B[0]) for x in range(len(A))]

            for i in range(len(A)):
                for j in range(len(B[0])):
                    for k in range(len(B)):
                        if (A[i][k] and B[k][j]):
                            for a in str(A[i][k]).split(' + '):
                                for b in str(B[k][j]).split(' + '):
                                    term = ('({0})*({1})').format(a,b)

                                    if res[i][j]:
                                        res[i][j] = res[i][j] + ' + {0}'.format(term)
                                    else:
                                        res[i][j] = term

            return res


        Qr_env = matrix_product(Lr, Qn)
        Qr_env = matrix_product(Qr_env, zip(*Lr))

        for i in range(N_ENV_STO):
            for j in range(N_ENV_STO):
                Qr[N_REAC+i][N_REAC+j] = Qr_env[i][j]

        #we fill Qr with Qc_dem and Qc_env
        for i in range(N_REAC):
            for j in range(N_REAC):
                Qr[i][j] = Qr_dem[i][j]

    
        #we split Ls into Ls_dem and Ls_env
        Ls_dem = [[0]*N_REAC for x in range(N_PAR_SV + N_OBS)]
        for i in range(N_PAR_SV + N_OBS):
            for j in range(N_REAC):
                Ls_dem[i][j] = Ls[i][j]

        Ls_env = [[0]*N_ENV_STO for x in range(N_PAR_SV + N_OBS)]
        for i in range(N_PAR_SV + N_OBS):
            for j in range(N_ENV_STO):
                Ls_env[i][j] = Ls[i][N_REAC + j]


        #####################################################################################
        ##we create 4 versions of Q (no_dem_sto, no_env_sto, no_dem_sto_no_env_sto and full)
        #####################################################################################

        Qs = matrix_product(Ls, Qr)
        Qs = matrix_product(Qs, zip(*Ls))

        Qs_dem = matrix_product(Ls_dem, Qr_dem)
        Qs_dem = matrix_product(Qs_dem, zip(*Ls_dem))

        Qs_env = matrix_product(Ls_env, Qr_env)
        Qs_env = matrix_product(Qs_env, zip(*Ls_env))

        calc_Q = {'no_dem_sto': {'Q_proc':[],
                                 'Q_obs':[],
                                 'Q': Qs_env},
                  'no_env_sto': {'Q_proc':[],
                                 'Q_obs':[],
                                 'Q': Qs_dem},
                  'full': {'Q_proc':[],
                           'Q_obs':[],
                           'Q': Qs},
                  'no_dem_sto_no_env_sto':{'Q_proc':[],
                                           'Q_obs':[],
                                           'Q': []}}

        if debug:
            for k in calc_Q:

                print '\n\nNon null term of Q_'+ k
                print "sv:"
                for i, x in enumerate(self.par_sv):
                    print i, x

                print "obs:"
                for i, x in enumerate(self.obs_var_def):
                    print N_PAR_SV+ i, x

                for i in range(len(calc_Q[k]['Q'])):
                    for j in range(i+1):
                        if calc_Q[k]['Q'][i][j]:
                            print '----------'
                            #print Q[i][j]
                            print 'Q[{0}][{1}]: '.format(i, j),  self.make_C_term(calc_Q[k]['Q'][i][j], True, human=True)
                            if i != j:
                                print 'Q[{0}][{1}] == Q[{1}][{0}]: '.format(i, j), self.make_C_term(calc_Q[k]['Q'][i][j], True, human=True) == self.make_C_term(calc_Q[k]['Q'][j][i], True, human=True)
        

        #convert in a version easy to template in C
        #Note that we only template the lower triangle (Q is symmetrical)
        for k, tpl in calc_Q.iteritems():
            if tpl['Q']:
                for i  in range(len(tpl['Q'])):
                    for j in range(i+1):
                        if tpl['Q'][i][j]:
                            if i< N_PAR_SV and j < N_PAR_SV:                        
                                tpl['Q_proc'].append({'i': i, 'j': j, 'rate': self.make_C_term(tpl['Q'][i][j], True)})
                            else:
                                tpl['Q_obs'].append({'i': {'is_obs': False, 'ind': i} if i < N_PAR_SV else {'is_obs': True, 'ind': i - N_PAR_SV},
                                                     'j': {'is_obs': False, 'ind': j} if j < N_PAR_SV else {'is_obs': True, 'ind': j - N_PAR_SV},
                                                     'rate': self.make_C_term(tpl['Q'][i][j], True)})

        ##cache special functions
        for key in calc_Q:            
            if calc_Q[key]['Q']:

                optim_rates_proc = [x['rate'] for x in calc_Q[key]['Q_proc']]
                optim_rates_obs = [x['rate'] for x in calc_Q[key]['Q_obs']]
                optim_rates = optim_rates_proc + optim_rates_obs

                calc_Q[key]['sf'] = self.cache_special_function_C(optim_rates, prefix='_sf[cac]')
                
                for i in range(len(optim_rates_proc)):
                    calc_Q[key]['Q_proc'][i]['rate'] = optim_rates[i]

                n_proc = len(optim_rates_proc)
                for i in range(len(optim_rates_obs)):
                    calc_Q[key]['Q_obs'][i]['rate'] = optim_rates[n_proc + i]

            else:
                calc_Q[key]['sf'] = []


        return calc_Q



if __name__=="__main__":

    """test Ccoder"""

    import json
    import os

    c = json.load(open(os.path.join('example', 'noise', 'context.json')))
    p = json.load(open(os.path.join('example', 'noise', 'process.json')))
    l = json.load(open(os.path.join('example', 'noise', 'link.json')))

    ##fix path (this is normally done by pmbuilder(1))
    for x in c['data'] + c['metadata']:
        x['source'] = os.path.join('example', 'noise', x['source'])
    

    model = Ccoder(c, p, l)
    
#    model.eval_Q(debug=True)

    caches = map(lambda x: model.make_C_term(x, True), ['sin(2*M_PI*(t/ONE_YEAR +r0))', 'sin(2*M_PI*(t/ONE_YEAR +r0))'])
    print model.cache_special_function_C(caches, prefix='_sf[cac]')
    print caches


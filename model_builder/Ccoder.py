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
from sympy import diff, Symbol, sympify
from sympy.printing import ccode
import copy

class Ccoder(Cmodel):
    """write the C code from the user input coming from the web interface..."""

    def __init__(self, context, process, link,  **kwargs):
        Cmodel.__init__(self, context, process, link,  **kwargs)

        self.set_par_fixed = set(self.par_fixed).union(set(['p_0', 'sum_SV', 'prop']))
        self.all_par = self.par_sv + self.drift_par_proc + self.par_proc +  self.par_obs + ['x'] + ['N']


    def toC(self, term, is_ode=False):

        if term in self.par_sv:
            return 'X[ORDER_{0}*N_CAC+cac]'.format(term)

        elif term in self.set_par_fixed:
            if term == 'p_0':
                return 'p_data->pop_size_t0[cac]'
            elif term == 'sum_SV':
                return  'sum_SV(X, cac)'
            elif term == 'prop':
                return 'p_data->rep1[nn][ts]'
            elif term in self.par_fixed_obs:
                return 'covar[ORDER_{0}][nn][ts]'.format(term)
            else:
                return 'covar[ORDER_{0}][nn][cac]'.format(term)

        elif term in self.par_proc :
            if term in self.drift_par_proc:
                return 'back_transform_x(X[p_data->drift[ORDER_drift__par_proc__{0}]->offset + routers[ORDER_{0}]->map[cac]], routers[ORDER_{0}]->map[cac], routers[ORDER_{0}])'.format(term)
            else:
                return 'par[ORDER_{0}][routers[ORDER_{0}]->map[cac]]'.format(term)

        elif term in self.par_obs:
            if term in self.drift_par_obs:
                return 'back_transform_x(X[p_data->drift[ORDER_drift__par_obs__{0}]->offset + routers[ORDER_{0}]->map[ts]], routers[ORDER_{0}]->map[ts], routers[ORDER_{0}])'.format(term)
            else:
                return 'par[ORDER_{0}][routers[ORDER_{0}]->map[ts]]'.format(term)

        elif term == 'correct_rate' and is_ode:
            return ''

        else: ##r is an operator or x
            return term


    def generator_C(self, term, is_ode=False):

        terms = self.change_user_input(term)

        ##if is_ode, remove operator after or before noise__ (noise__ term are not present in ODE models) (e.g *noise__ -> '' and noise__trans()* -> ''...
        if is_ode:
            noise = set([x for x in terms if 'noise__' in x])
            for n in noise:
                pos_starts = [i for i,x in enumerate(terms) if x == n]
                for pos_start in pos_starts:
                    pos_closing_bracket = pos_start
                    while terms[pos_closing_bracket] != ')':
                        pos_closing_bracket +=1

                    ##remove '*' operator after or before noise__xxx(xxx).
                    if (pos_start == 0) and terms[pos_closing_bracket+1] == '*': ##noise__ is the first term
                        terms[pos_closing_bracket+1] = ''

                    elif (pos_closing_bracket == (len(terms)-1)) and terms[pos_start-1] == '*': ##noise__ is the last term
                        terms[pos_start-1] = ''

                    elif (pos_start > 0) and (pos_closing_bracket != (len(terms)-1)): ##noise is surounded by other terms. In this case we remove the operator that is '*'
                        op_left = terms[pos_start-1]
                        op_right = terms[pos_closing_bracket+1]
                        if op_left == '*':
                            terms[pos_start-1] = ''
                        else:
                            terms[pos_closing_bracket+1] = ''

        ind = 0
        while (ind < len(terms)):
            if terms[ind].split('__')[0] in self.special_functions:
                myf = terms[ind].split('__')
                if myf[0] == 'noise':
                    if not is_ode:
                        yield self.toC(terms[ind], is_ode)
                    ##skip argument
                    while terms[ind] != ')':
                        ind +=1
                else :
                    while terms[ind] != ')':
                        yield self.toC(terms[ind], is_ode)
                        ind +=1

                ##add extra terms (no whitespace)
                if myf[0] == 'sinusoidal_forcing':
                    yield ',t'
                elif myf[0] == 'terms_forcing':
                    yield ',t,p_data,cac'
                elif myf[0] == 'step':
                    yield ',t'
                elif myf[0] == 'step_lin':
                    yield ',t'
                elif myf[0] == 'correct_rate' and not is_ode:
                    yield ',dt'

                ##close bracket
                if myf[0] != 'noise':
                    yield self.toC(terms[ind], is_ode)
                ind +=1

            else:
                yield self.toC(terms[ind], is_ode)
                ind += 1

    def make_C_term(self, term, is_ode=False, derivate=None):

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
            cterm = ccode(diff(sympify(safe), sy))
        else:
            cterm = ccode(sympify(safe))

        #remove the plom___ prefix
        term = cterm.replace('plom___', '')

        #make the plom C expression
        return ''.join(self.generator_C(term, is_ode=is_ode))


    def print_obs_prev(self):
        """generate C code to track the observed **prevalence** """

        Clist = []

        for i in range(len(self.obs_var_def)):
            if not isinstance(self.obs_var_def[i][0], dict): ##prevalence
                prev = ''
                for j in range(len(self.obs_var_def[i])):
                    prev += ' + p_X->proj[ORDER_{0}*N_CAC+c*N_AC+ac]'.format(self.obs_var_def[i][j])

                Clist.append(prev)

        return Clist


    def print_obs_inc_markov(self):
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
        for s in self.par_sv + self.universes:
            nbreac = len([r for r in self.proc_model if r['from']==s]) +1 ##+1 to stay in the same compartment or to declare smtg in case of no reaction (not super clean but makes C code easier...)
            Clist.append({'state':s, 'nb_reaction': nbreac})

        return Clist



    def get_gamma_noise_terms(self):
        gamma_noise_name = []
        gamma_noise_sd = []

        reac_noise = [self.change_user_input(n['rate']) for n in self.proc_model if 'noise__' in n['rate']]
        for terms in reac_noise:
            for x in terms:
                if 'noise__' in x and x not in gamma_noise_name:
                    gamma_noise_name.append(x)
                    ind = terms.index(x)
                    while terms[ind] != ')':
                        if terms[ind] in self.par_proc:
                            gamma_noise_sd.append(terms[ind])
                        ind +=1

        return zip(gamma_noise_name, gamma_noise_sd)


    def get_sd(self, rate):
        terms = self.change_user_input(rate)
        for x in terms:
            if 'noise__' in x:
                ind = terms.index(x)
                while terms[ind] != ')':
                    if terms[ind] in self.par_proc:
                        return terms[ind], x
                    ind+=1


    def make_prob(self, exitlist):
        probstring=''
        for r in exitlist:
            probstring+= ' + ' + r


    def print_prob(self):
        """general structure:

        sum=...;
        if(sum>LIKE_MIN){
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

        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        rates = list(set(r['rate'] for r in proc_model if r['from'] not in self.universes))
        caches = map(lambda x: self.make_C_term(x), rates)
        for r in proc_model:
            if r['from'] not in self.universes:
                r['ind_cache'] = rates.index(r['rate'])

        sf = self.cache_special_function_C(caches)

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

        return {'code': Ccode, 'caches': caches, 'sf': sf}


    def print_multinomial(self):
        Cstring=''
        for s in self.par_sv:
            nbexit = len([r for r in self.proc_model if r['from']==s])
            if nbexit>0:
                Cstring += 'sfr_ran_multinomial(p_calc->randgsl, {1}, (unsigned int) X[ORDER_{0}*N_CAC+cac], p_calc->prob[ORDER_{0}], p_calc->inc[ORDER_{0}][cac]);\n'.format(s,nbexit+1) ##+1 to stay in the compartment

        return Cstring


    def print_update(self):
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
                if myinput[nbreac]['to'] not in self.universes: ##we exclude deaths or transitions to DU in the update
                    incDict[myinput[nbreac]['to']] += ' + p_calc->inc[ORDER_{0}][cac][{1}]'.format(myinput[nbreac]['from'], nbreac)


        ##we add flow from the universes (Poisson term). We want to cache those flow so that the incidences can be computed
        poisson = []
        for s in self.universes:
            reac_from_univ = [r for r in self.proc_model if r['from'] == s]
            for nbreac in range(len(reac_from_univ)):
                poisson.append('p_calc->inc[ORDER_{0}][cac][{1}] = gsl_ran_poisson(p_calc->randgsl, ({2})*dt)'.format(s, nbreac, self.make_C_term(reac_from_univ[nbreac]['rate'])))
                incDict[reac_from_univ[nbreac]['to']] += ' + p_calc->inc[ORDER_{0}][cac][{1}]'.format(s, nbreac)


        Cstring=''
        for s in self.par_sv:
            Cstring += 'X[ORDER_{0}*N_CAC+cac] = {1};\n'.format(s, incDict[s])

        return {'poisson':poisson, 'Cstring': Cstring}


    def cache_special_function_C(self, caches_C, sf=None, prefix='_sf'):
        """caches_C: List of cached expression in C
        caches_C is modified in place
        sf: an optional list of unique special function to be cached
        returns sf (created if sf input is None)
        """

        if not sf:
            sf = []
            for term in caches_C:
                if any([x in term for x in self.cached]):
                    terms = self.change_user_input(term)
                    for x in terms:
                        if x in self.cached:
                            f = x + '('
                            ind = terms.index(x)+2 #skip first parenthesis
                            pos = 1 #counter for open parenthesis
                            while terms[ind] != ')' and pos > 0:
                                if terms[ind] == '(':
                                    pos += 1
                                    if terms[ind] == ')':
                                        pos -= 1

                                f += terms[ind]
                                ind +=1
                            f += terms[ind]

                            sf.append(f)

            sf = list(set(sf))

        for i, term in enumerate(caches_C):
            if any([x in term for x in self.cached]):
                for s in sf:
                    caches_C[i] = caches_C[i].replace(s, prefix + '[{0}]'.format(sf.index(s)))

        return sf


    def print_ode(self):

        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        odeDict = dict([(x, []) for x in self.par_sv])

        rates = list(set(r['rate'] for r in proc_model))
        caches = map(lambda x: self.make_C_term(x, True), rates)
        for i, r in enumerate(proc_model):
            r['ind_cache'] = rates.index(r['rate'])
            r['ind_dem_sto'] = i

        sf = self.cache_special_function_C(caches, prefix='_sf[cac]')

        def get_rhs_term(sign, cached, reaction):
            if 'noise__' in reaction['rate']:
                noise_sd, noise_name = self.get_sd(reaction['rate'])
                noise_sd = self.toC(noise_sd)
            else: 
                noise_name = None
                noise_sd= None
            
            return {'sign': sign, 'term': cached, 'noise_name': noise_name, 'noise_sd': noise_sd, 'ind_dem_sto': reaction['ind_dem_sto']}


        ################################
        ##Dynamic of the state variables
        ################################

        ##outputs
        for r in proc_model:
            if r['from'] not in self.universes:
                cached = '_r[cac][{0}]*X[ORDER_{1}{2}]'.format(r['ind_cache'], r['from'], '*N_CAC+cac')
                odeDict[r['from']].append(get_rhs_term('-', cached, r))

        ##inputs
        for r in proc_model:
            if r['to'] not in self.universes:
                if r['from'] not in self.universes:
                    cached = '_r[cac][{0}]*X[ORDER_{1}{2}]'.format(r['ind_cache'], r['from'], '*N_CAC+cac')
                else:
                    cached= '_r[cac][{0}]'.format(r['ind_cache'])

                odeDict[r['to']].append(get_rhs_term('+', cached, r))

        
        #######################################
        ##Dynamic of the observed **incidence**
        #######################################

        obs_list = []

        for i in range(len(self.obs_var_def)):
            true_ind_obs = str(i)
            eq = []

            if isinstance(self.obs_var_def[i][0], dict): ##incidence
                for j in range(len(self.obs_var_def[i])):
                    id_out = [proc_model.index(r) for r in proc_model if ((r['from'] == self.obs_var_def[i][j]['from']) and (r['to'] == self.obs_var_def[i][j]['to']) and (r['rate'] == self.obs_var_def[i][j]['rate'])) ]
                    for o in id_out:
                        reaction = proc_model[o]
                        if self.obs_var_def[i][j]['from'] in self.universes:
                            cached = '_r[cac][{0}]'.format(reaction['ind_cache'])
                        else:
                            cached = '_r[cac][{0}]*X[ORDER_{1}{2}]'.format(reaction['ind_cache'], self.obs_var_def[i][j]['from'], '*N_CAC+cac')

                        eq.append(get_rhs_term('+', cached, reaction))

                obs_list.append({'index':i, 'eq': eq})



        ##############################################################################################################
        ##we create the ODE and  4 versions of the SDE system (no_dem_sto, no_env_sto, no_dem_sto_no_env_sto and full)
        ##############################################################################################################
        unique_noises = self.get_gamma_noise_terms()
        unique_noises_names = [x[0] for x in unique_noises]

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

            #TODO get rid of the 'dt' for Euler Maruyama (should be handled on the C side as sqrt(dt))'
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

            #TODO get rid of the 'dt' for Euler Maruyama (should be handled on the C side as sqrt(dt))'
            func['ode']['obs'].append({'index': i, 'eq': eq})
            func['no_dem_sto_no_env_sto']['obs'].append({'index': i, 'eq': '({0})*dt'.format(eq)})
            func['no_dem_sto']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt {1}'.format(eq, env)})
            func['no_env_sto']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt {1}'.format(eq, dem)})
            func['full']['obs'].append({'index': myobs['index'], 'eq': '({0})*dt + {1} {2}'.format(eq, dem, env)})

        
        return {'func': func, 'caches': caches, 'sf': sf}




    def print_order(self):
        order_univ = []
        N_PAR_SV = len(self.par_sv)
        for i, X in enumerate(self.universes):
            order_univ.append({'name': X, 'order': N_PAR_SV+i})

        return {'var': self.par_sv + self.par_proc + self.par_obs, 'drift': self.drift_var, 'data': self.par_fixed, 'universe': order_univ}


    def print_like(self):

        ##WARNING right now only the discretized normal is supported.
        ##TODO: generalization for different distribution

        return {'mean':self.make_C_term(self.obs_model['mean']),
                'var':self.make_C_term(self.obs_model['var'])}



    def jac(self, sf_jac_only):
        """compute jacobian matrix of the process model (including
        observed variable) using Sympy


        sf_jac_only: list of cached special function generated by
        self.print_prob() used to get the index of caches_C for the
        jacobian matrix of simulation methods

        """

        my_model = copy.deepcopy(self.proc_model)
        odeDict = dict([(x,'') for x in self.par_sv])


        ##outputs
        for r in my_model:
            if r['from'] not in self.universes:
                rate= ' - (({0})*{1})'.format(r['rate'], r['from'])
                odeDict[r['from']] += rate

        ##inputs
        for r in my_model:
            if r['to'] not in self.universes:
                if r['from'] not in self.universes:
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
                    eq += odeDict[self.obs_var_def[i][j]]

            else: ##incidence
                for j in range(len(self.obs_var_def[i])):
                    id_out = [self.proc_model.index(r) for r in self.proc_model if ((r['from'] == self.obs_var_def[i][j]['from']) and (r['to'] == self.obs_var_def[i][j]['to']) and (r['rate'] == self.obs_var_def[i][j]['rate'])) ]
                    for o in id_out:
                        reaction = my_model[o]
                        if self.obs_var_def[i][j]['from'] in self.universes:
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
                Cterm = self.make_C_term(odeDict[self.par_sv[s]], is_ode=True, derivate=sy)
                jac[s].append(Cterm)
                jac_only[s].append(Cterm)
                caches.append(Cterm)
                caches_jac_only.append(Cterm)

            #see doc of kalman.c drift_derivative()
            for sy in self.drift_par_proc:
                Cterm = self.make_C_term(odeDict[self.par_sv[s]], is_ode=True, derivate=sy)
                jac_drift[s].append({'value': Cterm,
                                     'der': self.make_C_term(sy, is_ode=True),
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
                Cterm = self.make_C_term(obsList[o], is_ode=True, derivate=sy)
                jac_obs[o].append(Cterm)
                caches.append(Cterm)

            #see doc of kalman.c drift_derivative()
            for sy in self.drift_par_proc:
                Cterm = self.make_C_term(obsList[o], is_ode=True, derivate=sy)
                jac_obs_drift[o].append({'value': Cterm,
                                         'der': self.make_C_term(sy, is_ode=True),
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


    def jac_proc_obs(self):
        """compute jacobian matrix of the mean of the obs process (assumed to be Gaussian) using Sympy"""

        return self.make_C_term(self.obs_model['mean'], is_ode=True, derivate='x')



    def eval_Q(self):
        """create Ls and Qc
        
        Ls: Dispersion matrix of stochastic differential equation as
        defined is Sarkka 2006 phD.
        Ls is of size n*s with n == N_PAR_SV and s == number of noise terms
        
        Qc matrix of size s*s.

        Qc is composed of 2 blocks: Qc_dem (for demographic
        stochasticity) and Qc_sto (for environmental stochasticity)

        Qc_dem: diagonal matrix (one term per reaction: rate/N)

        Qc_env: we build it using the following procedure:

        - enumerate all noise terms coming from environmental sto.
        Let's say we have p noise terms, s reactions and n state
        variables.

        - build a "noise to rate" dispersion matrix Lc (s*p) that is
        analog to a soechiometric matrix (L above): each column
        illustrates the impact of a noise term on each reaction.  The
        component (i,j) of Lc is equal to the rate of reaction i if it
        is concerned with noise term j.

        - build a noise diffusion matrix Qc (p*p) that is diagonal,
        with diagonal terms equal to sto^2. 

        - Qc_env is then simply given by Lc Qn Lc'


        Note: we assume only one environmental noise term per reaction

        """
        proc_model = copy.deepcopy(self.proc_model) ##we are going to modify it...

        N_REAC = len(proc_model)
        N_PAR_SV = len(self.par_sv)
        N_OBS = len(self.obs_var_def)


        unique_noises = self.get_gamma_noise_terms()
        unique_noises_names = [x[0] for x in unique_noises]

        N_ENV_STO_UNIQUE = len(unique_noises)


        ##add sd and order properties to noisy reactions
        N_ENV_STO = 0
        for r in proc_model:
            if 'noise__' in r['rate']:
                r['sd'], noise_name = self.get_sd(r['rate'])
                r['order_env_sto_unique'] = unique_noises_names.index(noise_name)
                r['order_env_sto'] = N_ENV_STO
                N_ENV_STO += 1

        s = N_REAC + N_ENV_STO ##for demographic stochasticity, one independent noise term per reaction

        #we split Ls into Ls_proc and Ls_obs
        Ls_proc = [[0]*s for x in range(N_PAR_SV)]
        Ls_obs = [[0]*s for x in range(N_OBS)]

        diag_Qc_dem_tpl = []
        
        ###########################################
        # Create Ls_proc, Ls_obs and diag_Qc_dem  #
        ###########################################

        #state variables
        for B_dem_ind, r in enumerate(proc_model):
            is_noise = 'noise__' in r['rate']
            if is_noise:
                B_sto_ind = N_REAC + r['order_env_sto']

            if r['from'] not in self.universes:
                i = self.par_sv.index(r['from'])
                Ls_proc[i][B_dem_ind] -= 1 ##demographic stochasticity
                if is_noise:
                    Ls_proc[i][B_sto_ind] -= 1 ##env stochasticity

                Qc_term = '({0})*{1}'.format(r['rate'], r['from'])
            else:
                Qc_term = r['rate']

            if r['to'] not in self.universes:
                i = self.par_sv.index(r['to'])
                Ls_proc[i][B_dem_ind] += 1
                if is_noise:
                    Ls_proc[i][B_sto_ind] += 1

            diag_Qc_dem_tpl.append({'i': B_dem_ind, 'j': B_dem_ind, 'rate':self.make_C_term(Qc_term + '/' + self.myN, True)}) ##TODO: check if / self.myN is needed

        # observed variables 
        for i in range(len(self.obs_var_def)): #(for every obs variable)

            for B_dem_ind, r in enumerate(proc_model):
                is_noise = 'noise__' in r['rate']
                if is_noise:
                    B_sto_ind = N_REAC + r['order_env_sto']

                ##prevalence: (TODO: get rid of prevalence)
                if not isinstance(self.obs_var_def[i][0], dict):

                    # if it involves prevalence as input
                    if r['from'] in self.obs_var_def[i]:
                        Ls_obs[i][B_dem_ind] -= 1
                        if is_noise:
                            Ls_obs[i][B_sto_ind] -= 1

                    # if it involves prevalence as output
                    if r['to'] in self.obs_var_def[i]:
                        Ls_obs[i][B_dem_ind] += 1
                        if is_noise:
                            Ls_obs[i][B_sto_ind] += 1

                #incidence:
                else:
                    # for every incidence
                    for inc in self.obs_var_def[i]:
                        # if it involves incidence
                        if (r['from'] == inc['from']) and (r['to'] == inc['to']) and (r['rate'] == inc['rate']):
                            Ls_obs[i][B_dem_ind] += 1
                            if is_noise:
                                Ls_obs[i][B_sto_ind] += 1



        ############################
        ## Create Qc_env = Lc Qn Lc'
        ############################
        Lc = [[0]*N_ENV_STO_UNIQUE for x in range(N_ENV_STO)]
        Qn = [[0]*N_ENV_STO_UNIQUE for x in range(N_ENV_STO_UNIQUE)]

        for r in proc_model:
            if 'noise__' in r['rate']:
                if r['from'] not in self.universes:
                    Qn_term = '({0})*{1}'.format(r['rate'], r['from'])
                else:
                    Qn_term = r['rate']

                Lc[r['order_env_sto']][r['order_env_sto_unique']] = self.make_C_term(Qn_term, True) #True ensures that noise__ terms are removed from the rate (ODE)
                Qn[r['order_env_sto_unique']][r['order_env_sto_unique']] = 'pow({0}, 2)'.format(self.toC(r['sd']))

        
        #Qc_env = Lc Qn tLc:
        Qc_env = [[0]*N_ENV_STO for x in range(N_ENV_STO)]
    
        def matrix_product(A, B):
            res = [[0]*len(B[0]) for x in range(len(A))]

            for i in range(len(A)):
                for j in range(len(B[0])):
                    for k in range(len(B)):
                        if (A[i][k] and B[k][j]):                                        
                            term = ('({0})*({1})').format(A[i][k], B[k][j])

                            if res[i][j]: #should never happen
                               res[i][j] = res[i][j] + ' + {0}'.format(term)
                            else:
                               res[i][j] = term 

            return res


        Qc_env = matrix_product(Lc, Qn)
        Qc_env = matrix_product(Qc_env, zip(*Lc))

        #convert in a version easy to template in C
        Qc_env_tpl = []
        for i, row in enumerate(Qc_env):
            for j, term in enumerate(row):
                if Qc_env[i][j]:
                  Qc_env_tpl.append({'i': N_REAC + i, 'j': N_REAC + j, 'rate': term})  


        ################################################################################################
        ##we create 4 versions of Ls and Qc_diag (no_dem_sto, no_env_sto, no_dem_sto_no_env_sto and full
        ################################################################################################
        calc_Q = {'no_dem_sto': {'Ls_proc':[x[N_REAC:(N_REAC + N_ENV_STO)] for x in Ls_proc],
                                 'Ls_obs':[x[N_REAC:(N_REAC + N_ENV_STO)] for x in Ls_obs],
                                 'Qc_tpl': Qc_env_tpl},
                  'no_env_sto': {'Ls_proc':[x[0:N_REAC] for x in Ls_proc],
                                 'Ls_obs':[x[0:N_REAC] for x in Ls_obs],
                                 'Qc_tpl': diag_Qc_dem_tpl},
                  'full': {'Ls_proc': Ls_proc,
                           'Ls_obs': Ls_obs,
                           'Qc_tpl': copy.deepcopy(diag_Qc_dem_tpl + Qc_env_tpl)}, #we deepcopy as we are going to modify diag_Qc_dem_tpl and Qc_env_tpl
                  'no_dem_sto_no_env_sto': {'Ls_proc': [],
                                            'Ls_obs': [],
                                            'Qc_tpl': []}}

        #fix index for no_dem_sto
        for x in calc_Q['no_dem_sto']['Qc_tpl']:
            x['i'] -= N_REAC
            x['j'] -= N_REAC

        ##cache special functions
        for key in calc_Q:
            s = len(calc_Q[key]['Qc_tpl'])
            calc_Q[key]['s'] = s
            if s:
                optim_rates = [x['rate'] for x in calc_Q[key]['Qc_tpl']]
                calc_Q[key]['sf'] = self.cache_special_function_C(optim_rates, prefix='_sf')
                for i in range(len(optim_rates)):
                    calc_Q[key]['Qc_tpl'][i]['rate'] = optim_rates[i]
            else:
                calc_Q[key]['sf'] = []

        return calc_Q




if __name__=="__main__":

    """test Ccoder"""

    import json
    import os
    from Builder import PlomModelBuilder

    c = json.load(open(os.path.join('example', 'noise', 'context.json')))
    p = json.load(open(os.path.join('example', 'noise', 'process.json')))
    l = json.load(open(os.path.join('example', 'noise', 'link.json')))

    ##fix path (this is normally done by plom(1))
    for x in c['data']:
        x['source'] = os.path.join('example', 'noise', x['source'])

    model = PlomModelBuilder(os.path.join(os.getenv("HOME"), 'plom_test_model'), c, p, l)

    model.print_ode()

##    model.prepare()
##    model.write_settings()
##    model.code()
##    model.compile()

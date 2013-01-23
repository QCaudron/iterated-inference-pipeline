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
                return 'p_calc->natural_drifted_safe[ORDER_drift__par_proc__{0}][routers[ORDER_{0}]->map[cac]]'.format(term)
            else:
                return 'par[ORDER_{0}][routers[ORDER_{0}]->map[cac]]'.format(term)

        elif term in self.par_obs:
            if term in self.drift_par_obs:
                return 'p_calc->natural_drifted_safe[ORDER_drift__par_obs__{0}][routers[ORDER_{0}]->map[ts]]'.format(term)
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


    def print_build_markov(self):
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

        rates = list(set(r['rate'] for r in self.proc_model if r['from'] not in self.universes))
        caches = map(lambda x: self.make_C_term(x), rates)
        for r in self.proc_model:
            if r['from'] not in self.universes:
                r['ind_cache'] = rates.index(r['rate'])

        sf = self.cache_special_function_C(caches)

        Ccode=''

        for s in self.par_sv:
            myexit = [r for r in self.proc_model if r['from'] == s]
            exitlist=[]

            if len(myexit)>0:

                for e in myexit:
                    exitlist.append('_r[{0}]*DT'.format(e['ind_cache']))

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
                poisson.append('p_calc->inc[ORDER_{0}][cac][{1}] = gsl_ran_poisson(p_calc->randgsl, ({2})*DT)'.format(s, nbreac, self.make_C_term(reac_from_univ[nbreac]['rate'])))
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
        odeDict = dict([(x,'') for x in self.par_sv])

        rates = list(set(r['rate'] for r in self.proc_model))
        caches = map(lambda x: self.make_C_term(x, True), rates)
        for r in self.proc_model:
            r['ind_cache'] = rates.index(r['rate'])

        sf = self.cache_special_function_C(caches, prefix='_sf[cac]')

        ##outputs
        for r in self.proc_model:
            if r['from'] not in self.universes:
                rate= ' - _r[cac][{0}]*X[ORDER_{1}{2}]'.format(r['ind_cache'], r['from'], '*N_CAC+cac')
                odeDict[r['from']] += rate

        ##inputs
        for r in self.proc_model:
            if r['to'] not in self.universes:
                if r['from'] not in self.universes:
                    rate= ' + _r[cac][{0}]*X[ORDER_{1}{2}]'.format(r['ind_cache'], r['from'], '*N_CAC+cac')
                    odeDict[r['to']] += rate
                else:
                    rate= ' + _r[cac][{0}]'.format(r['ind_cache'])
                    odeDict[r['to']] += rate

        ##output the system...
        Cstring=''
        for s in self.par_sv:
            Cstring += 'f[ORDER_{0}{2}] = {1};\n'.format(s, odeDict[s], '*N_CAC+cac')

        #################
        ##generate C code to compute the dynamic of the observed **incidence** in case of ODE models and put in into a Clist
        ##Clist = [{'true_ind_obs':'i (i count both incidence and prevalence (N_OBS_ALL))', 'right_hand_side':'f(X[i,c,ac]'}]
        #################

        Clist = []

        for i in range(len(self.obs_var_def)):
            true_ind_obs = str(i)
            right_hand_side=''

            if isinstance(self.obs_var_def[i][0], dict): ##incidence

                for j in range(len(self.obs_var_def[i])):
                    id_out = [self.proc_model.index(r) for r in self.proc_model if ((r['from'] == self.obs_var_def[i][j]['from']) and (r['to'] == self.obs_var_def[i][j]['to']) and (r['rate'] == self.obs_var_def[i][j]['rate'])) ]
                    for o in id_out:
                        reaction = self.proc_model[o]
                        if self.obs_var_def[i][j]['from'] in self.universes:
                            right_hand_side += ' + _r[cac][{0}]'.format(reaction['ind_cache'])
                        else:
                            right_hand_side += ' + _r[cac][{0}]*X[ORDER_{1}{2}]'.format(reaction['ind_cache'], self.obs_var_def[i][j]['from'], '*N_CAC+cac')

                Clist.append({'true_ind_obs':true_ind_obs, 'right_hand_side':right_hand_side})

        return {'sys': Cstring, 'caches': caches, 'obs':Clist, 'sf': sf}


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
        """computes non-correlated noise term of the kalman Q matrix"""

        #results
        res = []

        #get reactions with noise
        for reac in self.proc_model:
            if 'noise__' in reac['rate']:
                sd = [] #list of noise intensities
                terms = self.change_user_input(reac['rate'])
                for x in terms:
                    if 'noise__' in x:
                        ind = terms.index(x)
                        while terms[ind] != ')':
                            if terms[ind] in self.par_proc:
                                sd.append(self.toC(terms[ind]))
                            ind +=1

                if reac['from'] not in self.universes:
                    rate = '({0})*{1}'.format(reac['rate'], reac['from'])
                else:
                    rate = reac['rate']

                res.append({'from': self.par_sv.index(reac['from']) if reac['from'] not in self.universes else self.par_sv.index(reac['to']),
                           'to': self.par_sv.index(reac['to']) if reac['to'] not in self.universes else self.par_sv.index(reac['from']),
                           'prod_sd': '*'.join(sd),
                           'rate' : self.make_C_term(rate, True)}) #note: True ensure that noise__ terms are removed from the rate

        rates = [x['rate'] for x in res]
        sf = self.cache_special_function_C(rates, prefix='_sf')

        for i, r in enumerate(res):
            res[i]['rate'] = rates[i]

        return {'Q':res, 'sf': sf}



    def __str__(self):
        tobeprinted = ''

        jacs = self.jac()
        for i in jacs:
            tobeprinted += '\n\n' + i + ':\n'
            tobeprinted += '\n'.join(['  |  '.join(line) for line in jacs[i]])

        jac_proc_obs = self.jac_proc_obs()
        tobeprinted += '\n\njac_proc_obs:\n'
        tobeprinted += jac_proc_obs

        return tobeprinted


    def stoichiometric(self):
        """compute the stoichiometric and force of infection matrices for demographic stochasticity"""

        # init stoichiometric matrix
        S_sv = [['0']*len(self.proc_model) for s in self.par_sv] # state var bloc
        S_ov = [['0']*len(self.proc_model) for o in self.obs_var_def] # obs var bloc
        # init force of infection matrix
        F = ['0']*len(self.proc_model)


        ###################
        # state variables #
        ###################

        # for every atomic reaction
        rInd=0
        for r in self.proc_model:
            # for every state variable
            svInd=0
            for sv in self.par_sv:
                # state variable in left hand => output
                if r['from'] == sv:
                    S_sv[svInd][rInd] = '-1'
                # state variable in right hand => input
                elif r['to'] == sv:
                    S_sv[svInd][rInd] = '+1'
                svInd += 1

            # fill force matrix
            if r['from'] not in self.universes:
                rate = '({0})*{1}'.format(r['rate'], r['from'])
            else:
                rate = r['rate']

            rate += '/' + self.myN

            F[rInd] = self.make_C_term(rate, True)

            rInd += 1


        ######################
        # observed variables #
        ######################

        # for every observed variable
        for oInd in range(len(self.obs_var_def)):

            ###################
            # prevalence case #
            ###################

            if not isinstance(self.obs_var_def[oInd][0], dict): ##prevalence:
                # for every reaction
                rInd=0
                for r in self.proc_model:

                    # if it involves prevalence as input
                    if r['from'] == self.obs_var_def[oInd]:
                        S_ov[oInd][rInd]= '+1'

                    # if it involves prevalence as output
                    elif r['to'] == self.obs_var_def[oInd]:
                        S_ov[oInd][rInd]= '-1'

                    rInd += 1

            ##################
            # incidence case #
            ##################

            else: #incidence
                # for every incidence
                for inc in self.obs_var_def[oInd]:
                    # for every reaction
                    rInd=0
                    for r in self.proc_model:

                        # if it involves incidence
                        if (r['from'] == inc['from']) and (r['to'] == inc['to']) and (r['rate'] == inc['rate']):
                            S_ov[oInd][rInd]= '+1'

                        rInd += 1


        return {'S_sv': S_sv, 'S_ov': S_ov, 'F': F, 'rnb':len(self.proc_model)}



if __name__=="__main__":

    """test PLoM model builder"""

    print('see Builder for test')

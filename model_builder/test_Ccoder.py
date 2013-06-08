from Ccoder import Ccoder
import unittest
import copy

import json
import os


class TestCcoder(unittest.TestCase):

    def setUp(self):

        c = json.load(open(os.path.join('example', 'drift', 'context.json')))
        p = json.load(open(os.path.join('example', 'drift', 'process.json')))
        l = json.load(open(os.path.join('example', 'drift', 'link.json')))

        ##fix path (this is normally done by pmbuilder(1))
        for x in c['data'] + c['metadata']:
            x['source'] = os.path.join('example', 'drift', x['source'])

        self.m_drift = Ccoder(c, p, l)


        c = json.load(open(os.path.join('example', 'noise', 'context.json')))
        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        l = json.load(open(os.path.join('example', 'noise', 'link.json')))

        ##fix path (this is normally done by pmbuilder(1))
        for x in c['data'] + c['metadata']:
            x['source'] = os.path.join('example', 'noise', x['source'])

        self.m_noise = Ccoder(c, p, l)

        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        p["white_noise"][0]["reaction"][0] = {'to':'S','from':'U'}
        self.m_noise2 = Ccoder(c, p, l)


        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        p["white_noise"][0]["reaction"][0] = {'to':'R','from':'I'}
        self.m_noise3 = Ccoder(c, p, l)

        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        p["white_noise"][0]["reaction"][0] = {'to':'U','from':'I'}
        self.m_noise4 = Ccoder(c, p, l)

        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        p["model"].append({"from": "R", "to": "U",  "rate": "mu_d", "comment":"death"})
        p["white_noise"][0]["reaction"][0] = {'to':'U','from':'R'}
        self.m_noise5 = Ccoder(c, p, l)

        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        p["white_noise"][0]["reaction"].append({'to':'R','from':'I'})
        self.m_noise6 = Ccoder(c, p, l)

        p = json.load(open(os.path.join('example', 'noise', 'process.json')))
        p["white_noise"].append({"reaction":[{'to':'R','from':'I'}],"sd":"sto"})
        self.m_noise7 = Ccoder(c, p, l)


    def test_make_C_term(self):
        terms = [
            {'x': 'mu_b*(1.0+v*sin((v/N+(mu_b)))) + r0', #input
             'h': 'mu_b*(v*sin(mu_b + v/N) + 1.0) + r0', #expected human output
             'c': 'covar[ORDER_mu_b][nn][cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]*sin(covar[ORDER_mu_b][nn][cac]+par[ORDER_v][routers[ORDER_v]->map[cac]]/covar[ORDER_N][nn][cac])+1.0)+drifted[ORDER_drift__par_proc__r0][cac]'}, #expected C output

            {'x': 'N-S-I+S+I',
             'h': 'N',
             'c': 'covar[ORDER_N][nn][cac]'},

            {'x': 'rep*(1.0-rep)*prop*x + (rep*phi*prop*x)**2',
             'h': 'pow(phi, 2)*pow(rep, 2)*pow(x, 2)*pow(prop, 2) + rep*x*prop*(-rep + 1.0)',
             'c': 'pow(par[ORDER_phi][routers[ORDER_phi]->map[ts]],2)*pow(par[ORDER_rep][routers[ORDER_rep]->map[ts]],2)*pow(x,2)*pow(covar[ORDER_prop][nn][ts],2)+par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*covar[ORDER_prop][nn][ts]*(-par[ORDER_rep][routers[ORDER_rep]->map[ts]]+1.0)'},
        ]
            
        for t in terms:
            self.assertEqual(self.m_drift.make_C_term(t['x'], False, human=True), t['h'])
            self.assertEqual(self.m_drift.make_C_term(t['x'], False, human=False), t['c'])


    def test_make_C_term_derivate(self):
        x = 'sin(2*M_PI*(t/ONE_YEAR +r0))'
        h = '2*M_PI*cos(2*M_PI*(r0 + t/ONE_YEAR))'
        c = '2*M_PI*cos(2*M_PI*(drifted[ORDER_drift__par_proc__r0][cac]+t/ONE_YEAR))'

        self.assertEqual(self.m_drift.make_C_term(x, False, human=True, derivate='r0'), h)
        self.assertEqual(self.m_drift.make_C_term(x, False, human=False, derivate='r0'), c)


    def test_make_C_term_skip_correct_rate(self):
        #correct_rate is only skipped for C code

        x = 'mu_b*(1.0+correct_rate(v)*sin((correct_rate(v)/N+(mu_b)))) + r0'
        c = 'covar[ORDER_mu_b][nn][cac]*((par[ORDER_v][routers[ORDER_v]->map[cac]])*sin(covar[ORDER_mu_b][nn][cac]+(par[ORDER_v][routers[ORDER_v]->map[cac]])/covar[ORDER_N][nn][cac])+1.0)+drifted[ORDER_drift__par_proc__r0][cac]'

        self.assertEqual(self.m_drift.make_C_term(x, True, human=False), c)


    def test_make_C_term_extra_terms(self):
        terms = [
            {'x': 'terms_forcing(v)', 
             'c': 'terms_forcing(par[ORDER_v][routers[ORDER_v]->map[cac]],t,p_data,cac)'}, 

            {'x': 'step(v, r0)', 
             'c': 'step(par[ORDER_v][routers[ORDER_v]->map[cac]],drifted[ORDER_drift__par_proc__r0][cac],t)'}, 

            {'x': 'step_lin(v, r0)', 
             'c': 'step_lin(par[ORDER_v][routers[ORDER_v]->map[cac]],drifted[ORDER_drift__par_proc__r0][cac],t)'}, 

            {'x': 'correct_rate(v)', 
             'c': 'correct_rate(par[ORDER_v][routers[ORDER_v]->map[cac]],dt)'}, 
        ]
            
        for t in terms:
            self.assertEqual(self.m_drift.make_C_term(t['x'], False, human=False), t['c'])


    def test_cache_special_function_C(self):

        caches = map(lambda x: self.m_drift.make_C_term(x, True), ['sin(2*M_PI*(t/ONE_YEAR +r0))', 'sin(2*M_PI*(t/ONE_YEAR +r0))', 'terms_forcing(v)*sin(2*M_PI*(t/ONE_YEAR +r0))*terms_forcing(v)'])
        sf = self.m_drift.cache_special_function_C(caches, prefix='_sf[cac]')

        self.assertEqual(sf, ['sin(2*M_PI*(drifted[ORDER_drift__par_proc__r0][cac]+t/ONE_YEAR))', 'terms_forcing(par[ORDER_v][routers[ORDER_v]->map[cac]],t,p_data,cac)'])
        self.assertEqual(caches, ['_sf[cac][0]', '_sf[cac][0]', 'pow(_sf[cac][1],2)*_sf[cac][0]'])
        
    def test_eval_Q(self):
        calc_Q = self.m_noise.eval_Q()

        # testing env sto only
        term = '((((r0/N*v*I)*S)*((sto)**2))*((r0/N*v*I)*S)))'
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][0],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][1],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][3],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][4],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][0],'((1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][1],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][3],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][4],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][0],'((1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][1],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][3],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][4],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][0],'((1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][1],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][3],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][4],'((1)*'+term+'*(1)')

        # testing dem sto only
        term1 = '(mu_b*N))'
        term2 = '((r0/N*v*I)*S))'
        term3 = '((correct_rate(v))*I))'
        term4 = '((mu_d)*S))'
        term5 = '((mu_d)*I))'
        self.assertEqual(calc_Q["no_env_sto"]["Q"][0][0],'((1)*'+term1+'*(1) + ((-1)*'+term2+'*(-1) + ((-1)*'+term4+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][0][1],'((-1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][0][2],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q"][0][3],'((-1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][0][4],'((-1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][1][0],'((1)*'+term2+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][1][1],'((1)*'+term2+'*(1) + ((-1)*'+term3+'*(-1) + ((-1)*'+term5+'*(-1)' )
        self.assertEqual(calc_Q["no_env_sto"]["Q"][1][2],'((-1)*'+term3+'*(1) + ((-1)*'+term5+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][1][3],'((1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][1][4],'((1)*'+term2+'*(1) + ((-1)*'+term3+'*(-1) + ((-1)*'+term5+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][2][0],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q"][2][1],'((1)*'+term3+'*(-1) + ((1)*'+term5+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][2][2],'((1)*'+term3+'*(1) + ((1)*'+term5+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][2][3],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q"][2][4],'((1)*'+term3+'*(-1) + ((1)*'+term5+'*(-1)')       
        self.assertEqual(calc_Q["no_env_sto"]["Q"][3][0],'((1)*'+term2+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][3][1],'((1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][3][2],0)
        self.assertEqual(calc_Q["no_env_sto"]["Q"][3][3],'((1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][3][4],'((1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][4][0],'((1)*'+term2+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][4][1],'((1)*'+term2+'*(1) + ((-1)*'+term3+'*(-1) + ((-1)*'+term5+'*(-1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][4][2],'((-1)*'+term3+'*(1) + ((-1)*'+term5+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][4][3],'((1)*'+term2+'*(1)')
        self.assertEqual(calc_Q["no_env_sto"]["Q"][4][4],'((1)*'+term2+'*(1) + ((-1)*'+term3+'*(-1) + ((-1)*'+term5+'*(-1)')

        #print(calc_Q["no_env_sto"]["Q_proc"])

    
    def test_eval_Q_tricky_cases(self):

        calc_Q = self.m_noise2.eval_Q()
        # testing env sto only for m_noise2 : WN on U->S
        term = '(((mu_b*N)*((sto)**2))*(mu_b*N)))'
        for i in range(5):
            for j in range(5):
                if i==0 and j == 0:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q"][i][j],'((1)*'+term+'*(1)')
                else:
                    self.assertEqual(calc_Q["no_dem_sto"]["Q"][i][j],0)


        calc_Q = self.m_noise3.eval_Q()
        # testing env sto only for m_noise3 : WN on I->R
        term = '((((correct_rate(v))*I)*((sto)**2))*((correct_rate(v))*I)))'
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][1],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][2],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][4],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][1],'((1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][2],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][4],'((1)*'+term+'*(-1)')       
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][1],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][2],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][4],'((-1)*'+term+'*(-1)')
        

        calc_Q = self.m_noise4.eval_Q()
        # testing env sto only for m_noise4 : WN on I->U
        term = '((((mu_d)*I)*((sto)**2))*((mu_d)*I)))'
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][1],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][2],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][4],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][1],'((1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][2],'((1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][4],'((1)*'+term+'*(-1)')       
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][1],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][4],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][1],'((-1)*'+term+'*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][2],'((-1)*'+term+'*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][4],'((-1)*'+term+'*(-1)')

        calc_Q = self.m_noise5.eval_Q()
        # testing env sto only for m_noise5 : WN on R->U
        for i in range(5):
            for j in range(5):
                self.assertEqual(calc_Q["no_dem_sto"]["Q"][i][j],0)

        calc_Q = self.m_noise6.eval_Q()
        # testing env sto only for m_noise6 : correlated WN on I->R and S->I
        term1 = '((r0/N*v*I)*S)'
        term2 = '((correct_rate(v))*I)'
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][0],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][1],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][2],'((-1)*(('+term1+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][3],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][4],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][0],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][1],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][2],'((1)*(('+term1+'*((sto)**2))*'+term2+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][3],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][4],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][0],'((1)*(('+term2+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][1],'((1)*(('+term2+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][2],'((1)*(('+term2+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][3],'((1)*(('+term2+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][4],'((1)*(('+term2+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')       
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][0],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][1],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][2],'((1)*(('+term1+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][3],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][4],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][0],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][1],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][2],'((1)*(('+term1+'*((sto)**2))*'+term2+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][3],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][4],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term1+'))*(1) + ((1)*(('+term1+'*((sto)**2))*'+term2+'))*(-1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')

        calc_Q = self.m_noise7.eval_Q()
        # testing env sto only for m_noise7 : uncorrelated WN on I->R and S->I
        term1 = '((r0/N*v*I)*S)'
        term2 = '((correct_rate(v))*I)'
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][0],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][1],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][3],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][0][4],'((-1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][0],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][1],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][2],'((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][3],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][1][4],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][0],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][1],'((1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][2],'((1)*(('+term2+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][3],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][2][4],'((1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')       
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][0],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][1],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][2],0)
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][3],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][3][4],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][0],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][1],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][2],'((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][3],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1)')
        self.assertEqual(calc_Q["no_dem_sto"]["Q"][4][4],'((1)*(('+term1+'*((sto)**2))*'+term1+'))*(1) + ((-1)*(('+term2+'*((sto)**2))*'+term2+'))*(-1)')


if __name__ == '__main__':
    unittest.main()

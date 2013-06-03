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

        self.m = Ccoder(c, p, l) 


    def test_make_C_term(self):
        terms = [
            {'x': 'mu_b*(1.0+v*sin((v/N+(mu_b)))) + r0', #input
             'h': 'mu_b*(v*sin(mu_b + v/N) + 1.0) + r0', #expected human output
             'c': 'covar[ORDER_mu_b][nn][cac]*(par[ORDER_v][routers[ORDER_v]->map[cac]]*sin(covar[ORDER_mu_b][nn][cac]+par[ORDER_v][routers[ORDER_v]->map[cac]]/covar[ORDER_N][nn][cac])+1.0)+drifted[ORDER_drift__par_proc__r0][cac]'}, #expected C output

            {'x': 'N-S-I+S+I',
             'h': 'N',
             'c': 'covar[ORDER_N][nn][cac]'},

            {'x': 'rep*(1.0-rep)*prop*x + (rep*phi*prop*x)**2',
             'h': 'rep*x*prop*(1.0*pow(phi, 2)*rep*x*prop - 1.0*rep + 1.0)',
             'c': 'par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*covar[ORDER_prop][nn][ts]*(1.0*pow(par[ORDER_phi][routers[ORDER_phi]->map[ts]],2)*par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*covar[ORDER_prop][nn][ts]-1.0*par[ORDER_rep][routers[ORDER_rep]->map[ts]]+1.0)'},
        ]
            
        for t in terms:
            self.assertEqual(self.m.make_C_term(t['x'], False, human=True), t['h'])
            self.assertEqual(self.m.make_C_term(t['x'], False, human=False), t['c'])


    def test_make_C_term_derivate(self):
        x = 'sin(2*M_PI*(t/ONE_YEAR +r0))'
        h = '2*M_PI*cos(2*M_PI*(r0 + t/ONE_YEAR))'
        c = '2*M_PI*cos(2*M_PI*(drifted[ORDER_drift__par_proc__r0][cac]+t/ONE_YEAR))'

        self.assertEqual(self.m.make_C_term(x, False, human=True, derivate='r0'), h)
        self.assertEqual(self.m.make_C_term(x, False, human=False, derivate='r0'), c)


    def test_make_C_term_skip_correct_rate(self):
        #correct_rate is only skipped for C code

        x = 'mu_b*(1.0+correct_rate(v)*sin((correct_rate(v)/N+(mu_b)))) + r0'
        c = 'covar[ORDER_mu_b][nn][cac]*((par[ORDER_v][routers[ORDER_v]->map[cac]])*sin(covar[ORDER_mu_b][nn][cac]+(par[ORDER_v][routers[ORDER_v]->map[cac]])/covar[ORDER_N][nn][cac])+1.0)+drifted[ORDER_drift__par_proc__r0][cac]'

        self.assertEqual(self.m.make_C_term(x, True, human=False), c)


    


if __name__ == '__main__':
    unittest.main()

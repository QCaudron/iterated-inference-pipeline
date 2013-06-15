from Cmodel import Cmodel
import unittest
import copy

class TestCmodel(unittest.TestCase):

    def setUp(self):
        m = {
            'state': [{'id':'S'}, {'id':'I', 'tag': ['remainder', 'infectious']}, {'id': 'R'}],
            'parameter': [{'id':'r0'}, {'id':'v'}, {'id':'l'}, {'id':'e'}, {'id':'d'}, {'id':'sto'}, {'id':'alpha'}, {'id':'mu_b'}, {'id':'mu_d'}, {'id':'vol'}, {'id':'g'}],
            'model': [ 
                {'from': 'U', 'to': 'S',  'rate': 'mu_b*N'},
                {'from': 'S', 'to': 'E',  'rate': 'r0/N*v*(1.0+e*sin_t(d))*I', 'tag': 'transmission'},
                
                {'from': 'E', 'to': 'I', 'rate': '(1-alpha)*correct_rate(l)'},
                ##Here we split the reaction from E->U as we only observe a subpart
                {'from': 'E', 'to': 'U',  'rate': 'alpha*correct_rate(l)'},
                {'from': 'E', 'to': 'U',  'rate': 'mu_d'},
                
                {'from': 'S', 'to': 'U',  'rate': 'mu_d'},
                {'from': 'R', 'to': 'S',  'rate': 'g'},
                {'from': 'I', 'to': 'R', 'rate': '(1-alpha)*correct_rate(v)'},
                {'from': 'I', 'to': 'U',  'rate': 'alpha*correct_rate(v) + mu_d'}
            ],
            'diffusion': [{'parameter':'r0', 'volatility': 'vol', 'drift': 0.0}],
            'white_noise': [{'reaction': [{'from':'S', 'to': 'E'}], 'sd': 'sto'}]
        }
        
        ##context elements needed for Cmodel
        c = {
            'data': {},
            'metadata': [{'id': 'mu_b'}, {'id': 'mu_d'}, {'id': 'N'}, {'id': 'prop'}]
        }

        ##link elements needed for Cmodel
        l = {
            'observed': [
                {'id': 'Prev', 'definition': ['I'], 'model_id': 'common'},
                {'id': 'SI', 'definition': ['S', 'I'], 'model_id': 'common'},
                ##we have to specify a rate to the incidence E->U as we only observed a subpart of this reaction
                {'id': 'Inc_out', 'definition': [{'from':'I', 'to':'R'}, {'from':'E', 'to':'U', 'rate': 'mu_d'}], 'model_id': 'common'},
                {'id': 'Inc_in', 'definition': [{'from':'S', 'to':'E'}], 'model_id': 'common'},
                {'id': 'Inc_weird', 'definition': [{'from':'R', 'to':'S'}], 'model_id': 'common'}
            ],
            'observation': [{'id': 'common', 
                             'parameter': [{'id': 'rep','comment': 'reporting rate'}, 
                                           {'id': 'phi',  'comment': 'over-dispertion'}],
                             'model': {'distribution': 'discretized_normal',
                                       'mean': 'rep*prop*x',
                                       'var': 'rep*(1.0-rep)*prop*x + (rep*phi*prop*x)**2'}}]
        }


        #other model with different remainder
        m2 = copy.deepcopy(m)
        m2['state'] = [{'id':'S'}, {'id':'I', 'tag': 'infectious'}, {'id': 'R', 'tag': 'remainder'}]

        #no remainder
        m3 = copy.deepcopy(m)
        m3['state'] = [{'id':'S'}, {'id':'I', 'tag': 'infectious'}, {'id': 'R'}]

        self.m = Cmodel(c, m, l) 
        self.m2 = Cmodel(c, m2, l) 
        self.m3 = Cmodel(c, m3, l) 


    def test_change_user_input(self):
        x = self.m.change_user_input('r0*2*correct_rate(v)')
        self.assertEqual(x, ['r0', '*', '2', '*', 'correct_rate', '(', 'v', ')'])


    def test_par_sv(self):
        self.assertEqual(set(self.m.par_sv), set(['S','R']))
        self.assertEqual(set(self.m2.par_sv), set(['S','I']))
        self.assertEqual(set(self.m3.par_sv), set(['S','I', 'R']))

    def test_remaider(self):
        self.assertEqual(self.m.remainder, 'I')
        self.assertEqual(self.m2.remainder, 'R')
        self.assertEqual(self.m3.remainder, None)

    def test_par_proc(self):
        self.assertEqual(set(self.m.par_proc), set(['r0', 'v', 'l', 'e', 'd', 'sto', 'alpha', 'vol', 'g']))

    def test_par_obs(self):
        self.assertEqual(set(self.m.par_obs), set(['rep','phi']))

    def test_par_fixed(self):
        self.assertEqual(set(self.m.par_fixed), set(['N', 'prop', 'mu_b', 'mu_d']))

    def test_par_fixed_obs(self):
        self.assertEqual(self.m.par_fixed_obs, set(['prop']))

    def test_drift_par_proc(self):
        self.assertEqual(set(self.m.drift_par_proc), set(['r0']))

    def test_vol_par_proc(self):
        self.assertEqual(set(self.m.vol_par_proc), set(['vol']))

    def test_drift_par_obs(self):
        self.assertEqual(self.m.drift_par_obs, [])

    def test_vol_par_obs(self):
        self.assertEqual(self.m.vol_par_obs, [])

    def test_drift_var(self):
        self.assertEqual(set(self.m.drift_var), set(['drift__par_proc__r0']))

    def test_proc_model(self):

        expected = [
            {'from': 'U', 'to': 'S',  'rate': 'mu_b*N'},
            {'from': 'S', 'to': 'E',  'rate': 'r0/N*v*(1.0+e*sin_t(d))*(N-S-R)', "tag": 'transmission', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}, #change  
            {'from': 'E', 'to': 'I', 'rate': '(1-alpha)*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'alpha*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'mu_d'},            
            {'from': 'S', 'to': 'U',  'rate': 'mu_d'},
            {'from': 'R', 'to': 'S',  'rate': 'g'},
            {'from': 'I', 'to': 'R', 'rate': '((1-alpha)*correct_rate(v))*(N-S-R)'}, #change
            {'from': 'I', 'to': 'U',  'rate': '(alpha*correct_rate(v)+mu_d)*(N-S-R)'} #change
        ]

        expected2 = [
            {'from': 'U', 'to': 'S',  'rate': 'mu_b*N'},
            {'from': 'S', 'to': 'E',  'rate': 'r0/N*v*(1.0+e*sin_t(d))*I', "tag": 'transmission', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}, #change  
            {'from': 'E', 'to': 'I', 'rate': '(1-alpha)*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'alpha*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'mu_d'},            
            {'from': 'S', 'to': 'U',  'rate': 'mu_d'},
            {'from': 'R', 'to': 'S',  'rate': '(g)*(N-S-I)'}, #change
            {'from': 'I', 'to': 'R', 'rate': '(1-alpha)*correct_rate(v)'},
            {'from': 'I', 'to': 'U',  'rate': 'alpha*correct_rate(v)+mu_d'}
        ]


        expected3 = [
            {'from': 'U', 'to': 'S',  'rate': 'mu_b*(S+I+R)'}, #change  
            {'from': 'S', 'to': 'E',  'rate': 'r0/(S+I+R)*v*(1.0+e*sin_t(d))*I', "tag": 'transmission', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}, #change  
            {'from': 'E', 'to': 'I', 'rate': '(1-alpha)*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'alpha*correct_rate(l)'},
            {'from': 'E', 'to': 'U',  'rate': 'mu_d'},            
            {'from': 'S', 'to': 'U',  'rate': 'mu_d'},
            {'from': 'R', 'to': 'S',  'rate': 'g'},
            {'from': 'I', 'to': 'R', 'rate': '(1-alpha)*correct_rate(v)'},
            {'from': 'I', 'to': 'U',  'rate': 'alpha*correct_rate(v)+mu_d'}
        ]



        self.assertEqual(self.m.proc_model, expected)
        self.assertEqual(self.m2.proc_model, expected2)
        self.assertEqual(self.m3.proc_model, expected3)


    def test_obs_var(self):
        self.assertEqual(self.m.obs_var, ['Inc_out', 'Inc_in', 'Inc_weird', 'Prev', 'SI'])

    def test_obs_var_def(self):
        expected = [
            [{'from': 'I', 'to': 'R', 'rate': '((1-alpha)*correct_rate(v))*(N-S-R)'}, {'from': 'E', 'to': 'U', 'rate': 'mu_d'}],
            [{'from': 'S', 'to': 'E', 'rate': 'r0/N*v*(1.0+e*sin_t(d))*(N-S-R)', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}],
            [{'from': 'R', 'to': 'S', 'rate': 'g'}],
            ['I'],
            ['S', 'I']
        ]

        expected2 = [
            [{'from': 'I', 'to': 'R', 'rate': '(1-alpha)*correct_rate(v)'}, {'from': 'E', 'to': 'U', 'rate': 'mu_d'}],
            [{'from': 'S', 'to': 'E', 'rate': 'r0/N*v*(1.0+e*sin_t(d))*I', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}],
            [{'from': 'R', 'to': 'S', 'rate': '(g)*(N-S-I)'}],
            ['I'],
            ['S', 'I']
        ]

        expected3 = [
            [{'from': 'I', 'to': 'R', 'rate': '(1-alpha)*correct_rate(v)'}, {'from': 'E', 'to': 'U', 'rate': 'mu_d'}],
            [{'from': 'S', 'to': 'E', 'rate': 'r0/(S+I+R)*v*(1.0+e*sin_t(d))*I', 'white_noise': {'name': 'white_noise__0', 'sd': 'sto'}}],
            [{'from': 'R', 'to': 'S', 'rate': 'g'}],
            ['I'],
            ['S', 'I']
        ]

        self.assertEqual(self.m.obs_var_def, expected)
        self.assertEqual(self.m2.obs_var_def, expected2)
        self.assertEqual(self.m3.obs_var_def, expected3)


if __name__ == '__main__':
    unittest.main()

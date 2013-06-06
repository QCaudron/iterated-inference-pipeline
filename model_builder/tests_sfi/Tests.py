import unittest
import os
import subprocess
import shutil
import json
from numpy import genfromtxt
import math

class TestLogTransfsAndPMCMC(unittest.TestCase):
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/noise/model')

            
      def test_prior_unif_transf_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_normal_transf_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            test = self.call_test_normal('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_normal_and_unif_transf_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96,r0:city2__all:prior:uniform,r0:city2__all:guess:10,r0:city2__all:min:9.5,r0:city2__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96,r0:city2__all:prior:uniform,r0:city2__all:guess:10,r0:city2__all:min:9.5,r0:city2__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            test1 = self.call_test_normal('r0.city1__all')
            test2 = self.call_test_unif('r0.city2__all')
            self.assertEqual(int(test1)*int(test2),1)

      def test_prior_unif_transf_logit_ab(self):
            os.system('plom pipe theta.json -S r0:transformation:logit_ab,r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:transformation:logit_ab,r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_unif_transf_identity(self):
            os.system('plom pipe theta.json -S r0:transformation:identity,r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:transformation:identity,r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_prior_normal_transf_identity(self):
            os.system('plom pipe theta.json -S r0:transformation:identity,r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:transformation:identity,r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            test = self.call_test_normal('r0.city1__all')
            self.assertEqual(test,'1')

      def call_test_unif(self,varname):
            shutil.copyfile(Root+'/TestsR/test_unif.r',Root+'/noise/model/test_R0_unif.r')
            os.system('R --vanilla < test_R0_unif.r ' + varname)
            f = open(Root+"/noise/model/outfile.txt","r")
            x = f.readlines()
            return x[0]

      def call_test_normal(self,varname):
            shutil.copyfile(Root+'/TestsR/test_normal.r',Root+'/noise/model/test_R0_normal.r')
            os.system('R --vanilla < test_R0_normal.r ' + varname)
            f = open(Root+"/noise/model/outfile.txt","r")
            x = f.readlines()
            return x[0]


class TestKalmanOnDiffusions(unittest.TestCase):
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/noise_test_diff/model')

      def test_1step(self):
            os.system('plom pipe theta.json | ./kalman -o 2')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['low95drifttest_parall'][0],-1.96,5)

      def test_10step(self):
            os.system('plom pipe theta.json | ./kalman -o 10')
            tab = genfromtxt('hat_0.csv',delimiter=',',names=True)
            self.assertAlmostEqual(tab['low95drifttest_parall'][9]/math.sqrt(10),-1.96,5)


def suite_TestLogTransfsAndPMCMC():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestLogTransfsAndPMCMC))
      return suite

def suite_TestKalmanOnDiffusions():
      suite = unittest.TestSuite()
      suite.addTest(unittest.makeSuite(TestKalmanOnDiffusions))
      return suite

if __name__ == '__main__' :

      run_LogTransfsAndPMCMC = 0
      run_KalmanOnDiffusions = 1
      

      Root = os.getcwd()

      if run_LogTransfsAndPMCMC:
            # copy noise from the examples and build it
            shutil.copytree(Root + '/../example/noise',Root + '/noise')     
            os.chdir(Root+'/noise')
            os.system('plom build -t theta.json --local')

            unittest.TextTestRunner().run(suite_TestLogTransfsAndPMCMC())

            shutil.rmtree(Root + '/noise')

      if run_KalmanOnDiffusions:
            # copy noise from the examples, add a diffusing parameter for tests, and build it
            shutil.copytree(Root + '/../example/noise',Root + '/noise_test_diff')     
            os.chdir(Root+'/noise_test_diff')
            p = json.load(open('process.json'))
            t = json.load(open('theta.json'))
            p["parameter"].append({'id':'test_par','comment':'will be freely diffusing'})
            p["parameter"].append({'id':'test_vol','comment':'will be equal to 1'})
            p["diffusion"]=[]
            p["diffusion"].append({'parameter':'test_par','volatility':'test_vol','drift':0.0})
            t["parameter"]['test_par']={'min':1,'max':1,'guess':1,'sd_transf':0.0}
            t["parameter"]['test_vol']={'min':1,'max':1,'guess':1,'sd_transf':0.0}
            json.dump(p,open('process.json','w'))
            json.dump(t,open('theta.json','w'))
            os.system('plom build -t theta.json --local')

            unittest.TextTestRunner().run(suite_TestKalmanOnDiffusions())

            shutil.rmtree(Root + '/noise_test_diff')


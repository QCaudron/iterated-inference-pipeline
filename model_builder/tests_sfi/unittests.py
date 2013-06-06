import unittest
import os
import subprocess
import shutil

class TestLogTransfsAndPMCMC(unittest.TestCase):
      def setUp(self):
      # Things that need to be done before tests
            os.chdir(Root+'/noise/model')

            
      def test_unif_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:city1__all:prior:uniform,r0:city1__all:guess:10,r0:city1__all:min:9.5,r0:city1__all:max:10.5 | ./pmcmc ode -M 10000 -S 1000 -E 1000 -o 0 --full')
            test = self.call_test_unif('r0.city1__all')
            self.assertEqual(test,'1')

      def test_normal_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            test = self.call_test_normal('r0.city1__all')
            self.assertEqual(test,'1')

      def test_normal_and_log_log(self):
            os.system('plom pipe theta.json -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96,r0:city2__all:prior:uniform,r0:city2__all:guess:10,r0:city2__all:min:9.5,r0:city2__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            os.system('plom pipe theta.json -B -C -S r0:city1__all:prior:normal,r0:city1__all:guess:10,r0:city1__all:min:8.04,r0:city1__all:max:11.96,r0:city2__all:prior:uniform,r0:city2__all:guess:10,r0:city2__all:min:9.5,r0:city2__all:max:10.5 | ./pmcmc ode -M 10000 -S 10000 -E 10000 -o 0 --full')
            test1 = self.call_test_normal('r0.city1__all')
            test2 = self.call_test_unif('r0.city2__all')
            self.assertEqual(int(test1)*int(test2),1)

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


if __name__ == '__main__' :

      Root = os.getcwd()
      ## Preliminary operations

      # copy noise from the examples and build it
      shutil.copytree(Root + '/../example/noise',Root + '/noise')     
      os.chdir(Root+'/noise')
      os.system('plom build -t theta.json --local')
      
      try:
            unittest.main()
      except:
            pass
      shutil.rmtree(Root + '/noise')

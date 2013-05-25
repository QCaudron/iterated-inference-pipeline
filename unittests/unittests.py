import unittest
import os
import subprocess

class TestTransfsAndPMCMC(unittest.TestCase):
      #def setUp(self):
      # Things that need to be done before tests

      def test_unif_log(self):
            Root = os.getcwd()
            print Root
            ## TODO: add command to recompile code
            os.chdir(Root+'/test_noise_unif_log/model')
            os.system('fit theta | ./pmcmc ode -M 1000 -o 0 --full')
            os.system('fit theta -B -C | ./pmcmc ode -M 1000 -o 0 --full')
            os.system('R --vanilla < R_test_unif.r')
            f = open(Root+'/test_noise_unif_log/model/' + "outfile.txt","r")
            x = f.readlines()
            os.chdir(Root)
            print(x[0])
            self.assertEqual(x[0],'1')

      def test_normal_log(self):
            Root = os.getcwd()
            print Root
            ## TODO: add command to recompile code
            os.chdir(Root+'/test_noise_normal_log/model')
            os.system('fit theta | ./pmcmc ode -M 1000 -o 0 --full')
            os.system('fit theta -B -C | ./pmcmc ode -M 1000 -o 0 --full')
            os.system('R --vanilla < R_test_normal.r')
            f = open(Root+'/test_noise_normal_log/model/' + "outfile.txt","r")
            x = f.readlines()
            os.chdir(Root)
            print(x[0])
            self.assertEqual(x[0],'1')
            

if __name__ == '__main__' :
      unittest.main()
             

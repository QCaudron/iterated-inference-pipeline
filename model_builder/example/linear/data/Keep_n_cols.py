

import csv
import os


name = 'prop.csv'
NbComumnsToKeep = 2

os.system('mv ' + name + ' tmp.csv')

with open('tmp.csv','rU') as csvfile:
    with open(name,'wb') as csvfile2:
    	 file = csv.reader(csvfile)
         file2 = csv.writer(csvfile2)

         ind = 0
	 for row in file:
             i=0
             for x in row[0:NbComumnsToKeep]:
                 if i > 0:
                      csvfile2.write(',')
                 if ind>0 and i == 1:
                     csvfile2.write(str(x))
                 else:
                     csvfile2.write('"'+str(x)+'"')
                 i=i+1
             csvfile2.write('\n')
             ind = ind+1

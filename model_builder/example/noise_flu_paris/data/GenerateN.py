

import csv


with open('data.csv','rU') as csvfile:
     with open('N.csv','wb') as csvfile2:


          file = csv.reader(csvfile)
          file2 = csv.writer(csvfile2,dialect='excel')

          ind = 0
          for row in file:
               if ind>0:
                    row[1] = 11700000


                    csvfile2.write('"'+row[0]+'"')
                    csvfile2.write(',')
                    csvfile2.write(str(row[1]))
                    csvfile2.write('\n')


               else:
                    row[1] = '"Paris__all"'

                    i = 0
                    for x in row:
                         if i > 0:
                              csvfile2.write(',')
                         csvfile2.write('"'+str(x)+'"')
                         i = i+1
                    csvfile2.write('\n')
          
               ind = ind + 1
	

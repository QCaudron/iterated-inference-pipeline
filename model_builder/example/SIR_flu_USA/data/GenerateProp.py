

import csv


with open('data.csv','rU') as csvfile:
     with open('prop.csv','wb') as csvfile2:


          file = csv.reader(csvfile)
          file2 = csv.writer(csvfile2,dialect='excel')

          ind = 0
          for row in file:
               if ind>0:
                    row[1] = 1
                    
                    csvfile2.write('"'+row[0]+'"')
                    csvfile2.write(',')
                    csvfile2.write(str(row[1]))
                    csvfile2.write('\n')

               else:
                    row[0] = '"date"'
                    row[1] = '"NJ__FluTrends__inc"'

                    i = 0
                    for x in row:
                         if i > 0:
                              csvfile2.write(',')
                         csvfile2.write('"'+str(x)+'"')
                         i = i+1
                    csvfile2.write('\n')
          
               ind = ind + 1

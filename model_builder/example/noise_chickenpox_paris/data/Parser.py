




import csv


with open('data_wrongformat.csv','rU') as csvfile:
    with open('data.csv','wb') as csvfile2:
    	 file = csv.reader(csvfile)
         file2 = csv.writer(csvfile2)

        
	 ind = 0
	 for row in file:
	     if ind>0:
                tmp = row[0].split("/")

                day = tmp[1]
                if len(day) == 1:
                   day = '0' + day

                month = tmp[0]
                if len(month) == 1:
                   month = '0' + month

                year = tmp[2]
                if year < 20:
                     year = '20'+year
                else:
                     year = '19'+year

                row[0] = '"' + year + "-" + month + "-" + day + '"'

                i = 0
                for x in row:
                   if i > 0:
                      csvfile2.write(',')
                   csvfile2.write(str(x))
                   i = i+1
                csvfile2.write('\n')
	     else:
                i = 0
                for x in row:
                   if i > 0:
                      csvfile2.write(',')
                   csvfile2.write(str(x))
                   i = i+1
                csvfile2.write('\n')

             ind = ind + 1

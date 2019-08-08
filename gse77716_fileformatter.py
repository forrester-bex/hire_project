# in bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77716/matrix/GSE77716_series_matrix.txt.gz


import pandas

# load file
inputfile = 'GSE77716_series_matrix.txt'

with open(inputfile, 'r') as inputf:
     counter = 1
     outfile = 'GSE42861_meth.txt'
     with open(outfile, 'w') as out:
             for line in inputf:
                     if counter >= 69 and counter <= 485648:
                             out.write(line+'\n')
                     counter += 1



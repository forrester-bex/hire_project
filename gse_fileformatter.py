inputfile = 'GSE42861_series_matrix.txt'

# read in file
with open(inputfile, 'r') as inputf:
	counter = 1
	# write out to file the lines associated with meth matrix
	outfile = 'GSE42861_meth.txt'
	with open(outfile, 'w') as out:
		for line in inputf:
			if counter >= 69 and counter <= 485648:
				out.write(line+'\n')
		
			counter += 1


meth_list = []

#formatting the data into a dataframe
import pandas
with open('GSE42861_meth.txt', 'r') as meth_data:
		for line in meth_data:
		# remove start and end lines
			if not line.startswith('!'):
				# split line by tabs
				splitlist = line.strip('\n').split('\t')
				# remove empty lines
				if len(splitlist) > 1:
				#append to new list in memory
					meth_list.append(splitlist)


# column names are row 1, ignoring first entry which references title of index column
column_names = meth_list[0][1:]
col_names = []
# strip of redundant ""
for entry in column_names:
	col_names.append(entry.strip("'").strip('"'))

# create list of row names (column 1 of every row)
row_names = []
r_names = []
for entry in meth_list[1:]:
	row_names.append(entry[0])

for entry in row_names:
	r_names.append(entry.strip("'").strip('"'))

# create matrix
matrix = []
for entry in meth_list[1:]:
	newent = []
	# add values as floats rather than strings
	for value in entry[1:]:
		newent.append(float(value))
	matrix.append(newent)


# create dataframe
df = pandas.DataFrame(matrix, columns=column_names, index=row_names)

# drop columns for samples GSM1051535 and GSM1051691
del df['GSM1051535']
del df['GSM1051691']

# write to csv file so can be easily loaded back in
df.to_csv('RAmeth_matrix.csv')


# calculate mean score of CpG loci (axis=1 means compute mean or row rather than column)
mean_scores = df.mean(axis=1)
mean_scores_list = []

for entry in mean_scores:
	mean_scores_list.append(entry)

ms = mean_scores.to_frame()
ms_csv = ms.to_csv()
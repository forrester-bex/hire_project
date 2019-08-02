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
with open('GSE42861_meth.txt', 'r') as meth_data:
		for line in meth_data:
		# remove start and end lines
			if not line.startswith('!'):
				# split line by tabs
				splitlist = line.strip('\n').split('\t')
				#append to new list in memory
				meth_list.append(splitlist)
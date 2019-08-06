inputfile = 'GSE42861_series_matrix.txt'


## CREATING THE METHYLATION MATRIX ## 
# read in file
with open(inputfile, 'r') as inputf:
     counter = 1
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
df = pandas.DataFrame(matrix, columns=col_names, index=r_names)

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

# create list of indexes to remove as <0.2 or > 0.8
ms = mean_scores.to_frame()
ms_csv = ms.to_csv()

pairs = []
test = ms_csv.split('\n')
for entry in test:
	if entry != '':
		splitentry = entry.split(',')
		if float(splitentry[1]) <0.2 or float(splitentry[1]) > 0.8:
			if len(splitentry[0]) > 3:
				try:
					pairs.append(splitentry[0])
				except ValueError:
					print entry

# remove lines and create new df 
newdf = df.drop(pairs)
newdf.to_csv('RAmeth_matrix_linesremoved.csv')



### CREATING THE SAMPLE SHEET ###

inputfile = 'GSE42861_series_matrix.txt'


# read in file
with open(inputfile, 'r') as inputf:
	counter = 1
	# write out to file the lines associated with the sample sheet
	outfile = 'GSE42861_samples.txt'
	with open(outfile, 'w') as out:
		for line in inputf:
			if counter >= 35 and counter < 69:
				out.write(line+'\n')
		
			counter += 1


attributes = []
with open('GSE42861_samples.txt', 'r') as samplesheet:
		for line in samplesheet:
				category = line.split('\n')
				cat = category[0]
				if len(cat) > 1:
					aslist = cat.split('\t')
					# create list of lists
					attributes.append(aslist)

### format attributes ###
ra_status = []
for entry in attributes[11]:
	if 'subject: Patient' in entry:
		ra_status.append(1)
	if 'subject: Normal' in entry:
		ra_status.append(0)

gender = []
for entry in attributes[13]:
	if 'gender: f' in entry:
		gender.append(0)
	if 'gender: m' in entry:
		gender.append(1)

smoking_history = []
for entry in attributes[14]:
	if 'smoking status: occasional' in entry:
		smoking_history.append((0,0,1))
	if 'smoking status: ex' in entry:
		smoking_history.append((1,0,0))
	if 'smoking status: never' in entry:
		smoking_history.append((0,0,0))
	if 'smoking status: na' in entry:
		smoking_history.append('na')
	if 'smoking status: current' in entry:
		smoking_history.append((0,1,0))

age = []
for entry in attributes[12]:
	if "age" in entry:
		newage = int(entry[6:].strip('"'))
		age.append(newage)


patient_id = []
for entry in attributes[1][1:]:
	new_id = entry.strip("'").strip('"')
	patient_id.append(new_id)


## also generate batch id ##
batch_id = []
for entry in attributes[31]:
	batch_id.append(entry.strip('"').split('/')[-1].split('_')[1])

column_names = patient_id 

row_names = ['ra status', 'age', 'gender', 'smoking']

matrix = [ra_status, age, gender, smoking_history]	

# create dataframe
samplesdf = pandas.DataFrame(matrix, columns=column_names, index=row_names)

# drop columns for samples GSM1051535 and GSM1051691
del samplesdf['GSM1051535']
del samplesdf['GSM1051691']

samplesdf.to_csv('samplesheet_colremoved.csv')		

# create batch_id table
batch_id_df = pandas.DataFrame(batch_id[1:], index=patient_id)
batch_id_df = batch_id_df.drop(['GSM1051535', 'GSM1051691'], axis=0)
df.to_csv('samplebatches.csv')

		

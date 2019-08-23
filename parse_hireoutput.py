import pandas 
import glob

df = pandas.read_csv("RA_HIREoutputmu.csv")

# add row names
matrix = pandas.read_csv('RAmeth_matrix_linesremoved_10Kmostvar.csv')


cpgs = matrix[matrix.columns[0]]

# add in cpg column
df['cpg'] = cpgs

# make column the row names
df = df.set_index("cpg")

# create separate dataframes for each cell type
celltypes = [df[df.columns[0]], df[df.columns[1]], df[df.columns[2]], df[df.columns[3]],
			 df[df.columns[4]], df[df.columns[5]]]


rankedlist = 'data_cleaned_v2.refactor.rankedlist.txt'
ranked = []

with open(rankedlist, 'r') as r:
	for line in r:
		cpg = line.strip('\n')
		ranked.append(cpg)

# create a list of the top 100 discriminatory cpgs
#top = ranked[:100]
top = ranked[:50]


# convert dfs to csv files
counter = 1
for entry in celltypes:
	entry.to_csv("celltype"+str(counter)+'.csv')
	counter += 1

# load csv files back in 

filelist = glob.glob('celltype[1-6].csv')

finaldata = []
for entry in filelist:
	datalines = []
	with open(entry, 'r') as data:
		for line in data:
			split_line = line.strip('\n').split(',')
			# if cpg in top 100 list, add to list
			if split_line[0] in top:
				datalines.append(split_line)
	finaldata.append(datalines)

# convert back to dfs and export to csvs
finalframes = []
for entry in finaldata:
	newdf = pandas.DataFrame(entry)
	newdf = newdf.set_index(0)
	finalframes.append(newdf)

counter = 1
for entry in finalframes:
	entry.to_csv("celltype"+str(counter)+'top100.csv')
	counter += 1

celltypes = []
cpgs = []
for entry in finaldata:
	data = []
	for cpg in entry:
		data.append(cpg[1])
		if cpg[0] not in cpgs:
			cpgs.append(cpg[0])
	celltypes.append(data)

# create joint matrix
joint_df = pandas.DataFrame(celltypes)
joint_df = joint_df.transpose()
joint_df['cpgs'] = cpgs
joint_df = joint_df.set_index('cpgs')
joint_df.columns=['celltype1', 'celltype2', 'celltype3', 'celltype4', 'celltype5', 'celltype6']
# export to csv
joint_df.to_csv("hire_allcelltypes_top100discrim.csv")
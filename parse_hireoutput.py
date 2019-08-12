import pandas 

df = pandas.read_csv("RA_HIRE_pvalues.csv")

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

sorted_dfs = []
for celltype in celltypes:
	sorted_celltype = celltype.sort_values(ascending=False)
	sorted_dfs.append(sorted_celltype[:100])

counter = 1
for entry in sorted_dfs:
	entry.to_csv("celltype"+str(counter)+'.csv')
	counter += 1
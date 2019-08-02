# HIRE requires a methylation matrix where one row corresponds to a CpG site
# and one column represents a sample

# geograbber
library(geograbi)
samples <- geograbi.get.samples("GSE35069")
meth <- geograbi.get.data("GSE35069")
data <- read.delim(myfile)


# read row names
head(rownames(data))

# read column names
head(colnames(data))
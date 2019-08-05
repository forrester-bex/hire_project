# Library processing

# geograbber (download files)
library(geograbi)
samples <- geograbi.get.samples("GSE35069")
meth <- geograbi.get.data("GSE35069")


# read row names
rownames(data)

# read column names
colnames(data)



### RA data set - GSE42861 ###
# Sample info:
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz

#meth data:
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz



### EDITING SAMPLE INFO ### 

### Removing samples GSM1051535 and GSM1051691 ###
# find row numbers of samples
which(samples["geo_accession"]=="GSM1051535") # returns col 11
which(samples["geo_accession"]=="GSM1051691") # returns col 167

# remove rows 11 and 167
samples2 <- samples[-c(11, 167),]


df$C <- ifelse(grepl("D", df$A), "yes", "no")


### RA DATASET INPUT FILES ###
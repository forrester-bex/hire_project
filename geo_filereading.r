# load in gzipped file
myfile <- gzfile('file.csv.gz', 'rt')

# read gzipped file

data <- read.csv(myfile, header=F)
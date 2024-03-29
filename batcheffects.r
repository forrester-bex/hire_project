# install sva 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sva")

# install test data
BiocManager::install("bladderbatch")

# load sva
library(sva)

# load in the data 
samples <- data.matrix(read.csv('samplesheet_colremoved.csv', row.names = 1))
meth_data <- data.matrix(read.csv('RAmeth_matrix_linesremoved.csv', row.names = 1))
batch_inf <- data.matrix(read.csv('samplebatches.csv', row.names = 1))

# transpose samples and turn from matrix into a dataframe
pheno <- data.frame(t(samples))

batch <- batch_inf[,1]

modcombat <- model.matrix(~1, data=pheno)

#returns matrix adjusted for batch
combat_edata <- ComBat(dat=meth_data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# write to csv
write.csv(combat_edata, "batcheffectsadjusted.csv")
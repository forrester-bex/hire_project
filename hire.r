#install HIRE
library("devtools")
install_github("XiangyuLuo/HIREewas")

# load HIRE
library("HIREewas")

# loads in samples as matrix with col1 as row names
samples <- data.matrix(read.csv('samplesheet_colremoved.csv', row.names = 1))
meth_data <- data.matrix(read.csv('RAmeth_matrix_linesremoved_10Kmostvar.csv', row.names = 1))


#launch HIRE
ret_list <- HIRE(meth_data, samples, num_celltype=6)

#a p-value matrix from the uniform distribution
#write output to pdf
    pvalues <- matrix(runif(600), 100, 6)
    #Visualize this p-value matrix
    pdf("RA_pvalues.pdf")
    riskCpGpattern(pvalues,
    main_title="An example", hc_row_ind = FALSE)
    dev.off()


  # write HIRE output to file
  capture.output(summary(ret_list), file = "RA_HIREoutput.txt")
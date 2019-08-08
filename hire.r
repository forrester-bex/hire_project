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


pdf("RA_casecontrolstatus.pdf")
riskCpGpattern(ret_list$pvalues[seq_len(100), c(2,1,3)],
main_title="Detected association pattern\n with disease status", hc_row_ind = FALSE) #c(2,1,3) was used because of the label switching
dev.off()
#age
#Visualize the association pattern with the age in the first 100 CpG sites
#write output to pdf
pdf("RA_age.pdf")
riskCpGpattern(ret_list$pvalues[seq_len(100), 6+c(2,1,3)],
    main_title="Detected association pattern\n with age", hc_row_ind = FALSE)
dev.off()

#a p-value matrix from the uniform distribution
#write output to pdf
    pvalues <- matrix(runif(600), 100, 6)
    #Visualize this p-value matrix
    pdf("RA_pvalues.pdf")
    riskCpGpattern(pvalues,
    main_title="An example", hc_row_ind = FALSE)
    dev.off()